#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

"""Parser module for ArrayExpress source."""

import logging
from typing import Dict, Any, Tuple, List, Optional
from sqlalchemy import select
from sqlalchemy.engine import URL
from ftplib import FTP

from ensembl.core.models import Gene as GeneORM

from ensembl.production.xrefs.parsers.BaseParser import BaseParser

class ArrayExpressParser(BaseParser):
    def run(self, args: Dict[str, Any]) -> Tuple[int, str]:
        source_id = args.get("source_id")
        species_id = args.get("species_id")
        species_name = args.get("species_name")
        xref_file = args.get("file", "")
        db_url = args.get("extra_db_url")
        ensembl_release = args.get("ensembl_release")
        xref_dbi = args.get("xref_dbi")
        verbose = args.get("verbose", False)

        if not source_id or not species_id:
            raise AttributeError("Missing required arguments: source_id and species_id")

        # Extract db connection parameters from file name
        project, db_user, db_host, db_port, db_name, db_pass = self.extract_params_from_string(
            xref_file, ["project", "user", "host", "port", "dbname", "pass"]
        )
        db_user = db_user or "ensro"
        db_port = db_port or "3306"

        # Get the species name(s)
        species_id_to_names = self.species_id_to_names(xref_dbi)
        if species_name:
            species_id_to_names.setdefault(species_id, []).append(species_name)

        if not species_id_to_names.get(species_id):
            return 0, "Skipped. Could not find species ID to name mapping"

        # Look up the species in ftp server and check if active
        species_lookup = self.get_species_from_ftp()
        if not self.is_arryaexpress_active(species_lookup, species_id_to_names[species_id], verbose):
            return 0, "Skipped. ArrayExpress source not active for species"

        species_name = species_id_to_names[species_id][0]

        # Connect to the appropriate arrayexpress db
        arrayexpress_db_url = self.get_arrayexpress_db_url(
            project, db_user, db_pass, db_host, db_port, db_name, species_name, ensembl_release, db_url, verbose
        )

        if not arrayexpress_db_url:
            raise AttributeError(
                "Could not find ArrayExpress DB. Missing or unsupported project value. Supported values: ensembl, ensemblgenomes."
            )

        xref_count = 0

        # Get data from arrayexpress db
        arrayexpress_data = self.get_arrayexpress_data(arrayexpress_db_url)

        # Add direct xref for every current gene found
        for row in arrayexpress_data:
            xref_id = self.add_xref(
                {
                    "accession": row.stable_id,
                    "label": row.stable_id,
                    "source_id": source_id,
                    "species_id": species_id,
                    "info_type": "DIRECT",
                },
                xref_dbi,
            )
            self.add_direct_xref(xref_id, row.stable_id, "gene", "", xref_dbi)
            xref_count += 1

        result_message = f"Added {xref_count} DIRECT xrefs"
        return 0, result_message

    def get_species_from_ftp(self) -> Dict[str, bool]:
        ftp_server = "ftp.ebi.ac.uk"
        ftp_dir = "pub/databases/microarray/data/atlas/bioentity_properties/ensembl"

        species_lookup = {}

        with FTP(ftp_server) as ftp:
            ftp.login("anonymous", "-anonymous@")
            ftp.cwd(ftp_dir)
            remote_files = ftp.nlst()

        for file in remote_files:
            species = file.split(".")[0]
            species_lookup[species] = True

        return species_lookup

    def is_arryaexpress_active(self, species_lookup: Dict[str, bool], names: List[str], verbose: bool) -> bool:
        for name in names:
            if species_lookup.get(name):
                if verbose:
                    logging.info(f"Found ArrayExpress has declared the name {name}. This was an alias")
                return True
        return False

    def get_arrayexpress_db_url(self, project: str, db_user: str, db_pass: str, db_host: str, db_port: str, db_name: str, species_name: str, ensembl_release: str, db_url: str, verbose: bool) -> Optional[URL]:
        if db_host:
            return URL.create("mysql", db_user, db_pass, db_host, db_port, db_name)
        elif project == "ensembl":
            if verbose:
                logging.info("Looking for db in mysql-ens-sta-1")
            registry = "ensro@mysql-ens-sta-1:4519"
            return self.get_db_from_registry(species_name, "core", ensembl_release, registry)
        elif project == "ensemblgenomes":
            if verbose:
                logging.info("Looking for db in mysql-eg-staging-1 and mysql-eg-staging-2")
            registry = "ensro@mysql-eg-staging-1.ebi.ac.uk:4160"
            sta_db_url = self.get_db_from_registry(species_name, "core", ensembl_release, registry)
            if not sta_db_url:
                registry = "ensro@mysql-eg-staging-2.ebi.ac.uk:4275"
                return self.get_db_from_registry(species_name, "core", ensembl_release, registry)
            return sta_db_url
        elif db_url:
            return db_url

        return None

    def get_arrayexpress_data(self, arrayexpress_db_url: URL) -> List[Dict[str, Any]]:
        db_engine = self.get_db_engine(arrayexpress_db_url)
        with db_engine.connect() as arrayexpress_dbi:
            query = select(GeneORM.stable_id).where(
                GeneORM.biotype != "LRG_gene", GeneORM.is_current == 1
            )
            result = arrayexpress_dbi.execute(query).mappings().all()
        
        return result
