#  **********************************************************************
#  Copyright (C) 2020 Johns Hopkins University Applied Physics Laboratory
#
#  All Rights Reserved.
#  For any other permission, please contact the Legal Office at JHU/APL.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#  **********************************************************************

import ast
import json
import os
import shlex
import shutil
import zipfile
from datetime import datetime

import pydash
import sh
from bson import ObjectId
from flask import Blueprint, request
from pymodm import connect

from shared.config import config
from shared.log import logger
from system.controllers import classification_job, evaluation_job, job_queue, simulation_job, user, user_job
from system.extensions import FlaskExtensions
from system.job_queue_manager import push_job, set_job_manager_running_status, restart_job_queue_watchdog
from system.models.job_manager import JobType, JobMode
from system.utils.biocontainers import get_biocontainers, kill_running_container
from system.utils.encoder import json_encoder
from system.utils.readtypes import get_read_types
from system.utils.security import allowed_file

logger.info("{}".format(config.MONGO_URI))
connect(config.MONGO_URI)

jobs_bp = Blueprint("jobs", __name__, url_prefix=config.SERVER_API_CHROOT)

mongodb = FlaskExtensions.mongodb


def validate_files(zipped_paths):
    validation_checks = list()
    for path in zipped_paths:
        if os.path.splitext(path)[1] == ".tsv":
            validation_checks.append(os.path.basename(path) + " is " + tsv_validation(path))
        else:
            validation_checks.append(os.path.basename(path) + " is " + fastq_validation(path))
    return validation_checks


def tsv_validation(tsv_path):
    # Check if abundance profile sums to 1.00
    FIRST_VALIDATION = "'{{sum+=$2}} END{{printf sum}}' {}".format(tsv_path)
    SECOND_VALIDATION = "'{{print NF}}' {}".format(tsv_path)
    THIRD_VALIDATION = "'NF!=3 {{print NR}}' {}".format(tsv_path)
    sort = "-nu"
    head = "-n 1"
    v1 = float(sh.awk(shlex.split(FIRST_VALIDATION)))
    v2 = float(sh.head(sh.sort(sh.awk(shlex.split(SECOND_VALIDATION)), sort), head))
    v3 = sh.awk(shlex.split(THIRD_VALIDATION))

    if v3 != '' or v2 != 3:
        logger.warning("NOT IN TSV FORMAT: FILE MUST HAVE EXACTLY THREE COLUMNS. ROW {} IS MISSING "
                       "FIELDS".format(v3))
        return "NOT IN TSV FORMAT: FILE MUST HAVE EXACTLY THREE COLUMNS. ROW {} IS MISSING FIELDS".format(v3)

    logger.info("ABUNDANCE PROFILE SUMS TO {}".format(v1))
    if v1 != 1:
        logger.warning("NOT IN TSV FORMAT: 2ND COLUMN MUST SUM TO 1")
        return "NOT IN TSV FORMAT: 2ND COLUMN MUST SUM TO 1"

    logger.info("TSV IS IN PROPER FORMAT")
    return "Valid"


def fastq_validation(fastq_path):
    # ONLY looking at first 4 lines
    FIRST_VALIDATION = "NR==1 {}".format(fastq_path)
    SECOND_VALIDATION = "NR==3 {}".format(fastq_path)
    CUT = "-c1"
    THIRD_VALIDATION = "NR==2 {}".format(fastq_path)
    allowed = "[A\n|G\n|T\n|C\n|N\n]"

    v1 = str(sh.cut(sh.awk(shlex.split(FIRST_VALIDATION)), CUT))
    v2 = str(sh.cut(sh.awk(shlex.split(SECOND_VALIDATION)), CUT))
    v3 = all(char in allowed for char in str(sh.awk(shlex.split(THIRD_VALIDATION))))

    if (v1 != '@\n'):
        logger.info("NOT IN FASTQ FORMAT: FILE MUST START WITH '@'")
        return "NOT IN FASTQ FORMAT: FILE MUST START WITH '@'"
    if (v2 != '+\n'):
        logger.info("NOT IN FASTQ FORMAT: Line 4 MUST START WITH '+'")
        return "NOT IN FASTQ FORMAT: Line 4 MUST START WITH '+'"
    if (v3 is False):
        logger.info("NOT IN FASTQ FORMAT: LINE 2 MUST ONLY CONTAIN A,C,G,T, OR N")
        return "NOT IN FASTQ FORMAT: LINE 2 MUST ONLY CONTAIN A,C,G,T, OR N"

    logger.info("FASTQ IS IN PROPER FORMAT")
    return "Valid"


@jobs_bp.route("/jobs/all", methods=["GET"])  # see all pending jobs
def get_jobs():
    res = user_job.find_all(as_json=True)
    return json.dumps(res), 200


@jobs_bp.route("/jobs/unhidden", methods=["GET"])  # see unhidden pending jobs
def get_unhidden_jobs():
    res = user_job.find_unhidden_jobs(as_json=True)
    return json.dumps(res), 200


@jobs_bp.route("/jobs/hidden", methods=["GET"])  # see hidden pending jobs
def get_hidden_jobs():
    res = user_job.find_hidden_jobs(as_json=True)
    return json.dumps(res), 200
