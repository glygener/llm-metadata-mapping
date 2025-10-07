from __future__ import annotations

import argparse
import csv
import sys
import time
from typing import List, Tuple, Optional

from Bio import Entrez

def init_entrez(email: str, api_key: Optional[str]) -> None:
    """

    :type api_key: Optional[str]
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
