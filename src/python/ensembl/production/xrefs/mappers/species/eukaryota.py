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

"""Mapper extension module for species eukaryota."""

import logging
from typing import Dict, List, Tuple
from sqlalchemy.orm import aliased
from sqlalchemy import select, update, func, delete
from sqlalchemy.sql.expression import Select
from sqlalchemy.dialects.mysql import insert

from ensembl.xrefs.xref_update_db_model import (
    Source as SourceUORM,
    Xref as XrefUORM,
    DependentXref as DependentXrefUORM,
    ObjectXref as ObjectXrefUORM
)

from ensembl.core.models import (
    Gene as GeneORM,
    Transcript as TranscriptORM,
    Xref as XrefCORM,
    ExternalDb as ExternalDbORM,
    ObjectXref as ObjectXrefCORM
)

from ensembl.production.xrefs.mappers.BasicMapper import BasicMapper

class eukaryota(BasicMapper):
    def gene_display_xref_sources(self) -> Tuple[List[str], Dict[str, Select]]:
        sources_list = [
            "TAIR_SYMBOL",
            "RFAM",
            "RNAMMER",
            "TRNASCAN_SE",
            "Uniprot_gn",
            "ENA_GENE",
            "BROAD_U_maydis",
            "BROAD_F_oxysporum",
            "BROAD_G_zeae",
            "BROAD_G_moniliformis",
            "BROAD_P_infestans",
            "phyra_jgi_v1.1",
            "physo1_jgi_v1.1",
            "phatr_jgi_v2",
            "phatr_jgi_v2_bd",
            "PGD_GENE",
            "Mycgr3_jgi_v2.0_gene",
            "BROAD_Magnaporthe_DB",
            "PHYTOZOME_GMAX_GENE",
        ]

        ignore_queries = {}

        # Ignore EntrezGene labels dependent on predicted RefSeqs
        MasterXref = aliased(XrefUORM)
        DependentXref = aliased(XrefUORM)
        MasterSource = aliased(SourceUORM)
        DependentSource = aliased(SourceUORM)

        query = select(ObjectXrefUORM.object_xref_id.distinct()).where(
            ObjectXrefUORM.xref_id == DependentXrefUORM.dependent_xref_id,
            ObjectXrefUORM.master_xref_id == DependentXrefUORM.master_xref_id,
            DependentXrefUORM.dependent_xref_id == DependentXref.xref_id,
            DependentXrefUORM.master_xref_id == MasterXref.xref_id,
            MasterXref.source_id == MasterSource.source_id,
            DependentXref.source_id == DependentSource.source_id,
            MasterSource.name.like("Refseq%predicted"),
            DependentSource.name.like("EntrezGene"),
            ObjectXrefUORM.ox_status == "DUMP_OUT",
        )
        ignore_queries["EntrezGene"] = query

        query = (
            select(ObjectXrefUORM.object_xref_id)
            .join(XrefUORM, XrefUORM.xref_id == ObjectXrefUORM.xref_id)
            .join(SourceUORM, SourceUORM.source_id == XrefUORM.source_id)
            .where(
                ObjectXrefUORM.ox_status == "DUMP_OUT",
                XrefUORM.label.regexp_match("^LOC[[:digit:]]+"),
            )
        )
        ignore_queries["LOC_prefix"] = query

        return sources_list, ignore_queries

    def transcript_display_xref_sources(self) -> Tuple[List[str], Dict[str, Select]]:
        sources_list = [
            "RFAM",
            "RNAMMER",
            "TRNASCAN_SE",
            "Uniprot_gn_trans_name",
            "ENA_GENE",
            "BROAD_U_maydis",
            "BROAD_F_oxysporum",
            "BROAD_G_zeae",
            "BROAD_G_moniliformis",
            "BROAD_P_infestans",
            "phyra_jgi_v1.1",
            "physo1_jgi_v1.1",
            "phatr_jgi_v2",
            "phatr_jgi_v2_bd",
            "PGD_GENE",
            "Mycgr3_jgi_v2.0_gene",
            "BROAD_Magnaporthe_DB",
            "PHYTOZOME_GMAX_GENE",
        ]

        ignore_queries = {}

        # Ignore EntrezGene labels dependent on predicted RefSeqs
        MasterXref = aliased(XrefUORM)
        DependentXref = aliased(XrefUORM)
        MasterSource = aliased(SourceUORM)
        DependentSource = aliased(SourceUORM)

        query = select(ObjectXrefUORM.object_xref_id.distinct()).where(
            ObjectXrefUORM.xref_id == DependentXrefUORM.dependent_xref_id,
            ObjectXrefUORM.master_xref_id == DependentXrefUORM.master_xref_id,
            DependentXrefUORM.dependent_xref_id == DependentXref.xref_id,
            DependentXrefUORM.master_xref_id == MasterXref.xref_id,
            MasterXref.source_id == MasterSource.source_id,
            DependentXref.source_id == DependentSource.source_id,
            MasterSource.name.like("Refseq%predicted"),
            DependentSource.name.like("EntrezGene"),
            ObjectXrefUORM.ox_status == "DUMP_OUT",
        )
        ignore_queries["EntrezGene"] = query

        query = (
            select(ObjectXrefUORM.object_xref_id)
            .join(XrefUORM, XrefUORM.xref_id == ObjectXrefUORM.xref_id)
            .join(SourceUORM, SourceUORM.source_id == XrefUORM.source_id)
            .where(
                ObjectXrefUORM.ox_status == "DUMP_OUT",
                XrefUORM.label.regexp_match("^LOC[[:digit:]]+"),
            )
        )
        ignore_queries["LOC_prefix"] = query

        return sources_list, ignore_queries

    def gene_description_sources(self) -> List[str]:
        sources_list = [
            "TAIR_LOCUS",
            "PomBase_GENE",
            "PomBase_TRANSCRIPT",
            "Uniprot/SWISSPROT",
            "Uniprot/SPTREMBL",
            "BROAD_U_maydis",
            "BROAD_F_oxysporum",
            "BROAD_G_zeae",
            "BROAD_G_moniliformis",
            "BROAD_P_infestans",
            "phyra_jgi_v1.1",
            "physo1_jgi_v1.1",
            "phatr_jgi_v2",
            "phatr_jgi_v2_bd",
            "PGD_GENE",
            "BROAD_Magnaporthe_DB",
            "PGSC_GENE",
            "PHYTOZOME_GMAX_GENE",
            "RFAM",
            "TRNASCAN_SE",
            "RNAMMER",
        ]

        return sources_list

    def set_transcript_names(self) -> None:
        logging.info("Assigning transcript names from gene names")

        core_dbi = self.core().connect()

        # Reset transcript display xrefs
        core_dbi.execute(update(TranscriptORM).values(display_xref_id=None))

        # Get the max xref and object_xref IDs
        xref_id = core_dbi.execute(select(func.max(XrefCORM.xref_id))).scalar()
        xref_id = int(xref_id)
        object_xref_id = core_dbi.execute(
            select(func.max(ObjectXrefCORM.object_xref_id))
        ).scalar()
        object_xref_id = int(object_xref_id)

        # Get all genes with set display_xref_id
        query = select(
            GeneORM.gene_id,
            ExternalDbORM.db_name,
            XrefCORM.dbprimary_acc,
            XrefCORM.display_label,
            XrefCORM.description,
        ).where(
            GeneORM.display_xref_id == XrefCORM.xref_id,
            XrefCORM.external_db_id == ExternalDbORM.external_db_id,
        )
        for row in core_dbi.execute(query).mappings().all():
            # Get the ID of transcript name external DB
            external_db_id = core_dbi.execute(
                select(ExternalDbORM.external_db_id).where(
                    ExternalDbORM.db_name.like(f"{row.db_name}_trans_name")
                )
            ).scalar()

            if not external_db_id:
                raise LookupError(
                    f"No external_db_id found for '{row.db_name}_trans_name'"
                )

            # Get transcripts related to current gene
            query = (
                select(TranscriptORM.transcript_id)
                .where(TranscriptORM.gene_id == row.gene_id)
                .order_by(TranscriptORM.seq_region_start, TranscriptORM.seq_region_end)
            )
            for transcript_row in core_dbi.execute(query).mappings().all():
                object_xref_id += 1

                # Check if xref already exists
                insert_xref_id = core_dbi.execute(
                    select(XrefCORM.xref_id).where(
                        XrefCORM.external_db_id == external_db_id,
                        XrefCORM.display_label == row.display_label,
                        XrefCORM.version == 0,
                        XrefCORM.description == row.description,
                        XrefCORM.info_type == "MISC",
                        XrefCORM.info_text == "via gene name",
                    )
                ).scalar()

                if not insert_xref_id:
                    xref_id += 1

                    # Insert new xref
                    core_dbi.execute(
                        insert(XrefCORM)
                        .values(
                            xref_id=xref_id,
                            external_db_id=external_db_id,
                            dbprimary_acc=row.display_label,
                            display_label=row.display_label,
                            version=0,
                            description=row.description,
                            info_type="MISC",
                            info_text="via gene name",
                        )
                        .prefix_with("IGNORE")
                    )

                    insert_xref_id = xref_id

                # Insert object xref
                core_dbi.execute(
                    insert(ObjectXrefCORM).values(
                        object_xref_id=object_xref_id,
                        ensembl_id=transcript_row.transcript_id,
                        ensembl_object_type="Transcript",
                        xref_id=insert_xref_id,
                    )
                )

                # Set transcript dispay xref
                core_dbi.execute(
                    update(TranscriptORM)
                    .values(display_xref_id=insert_xref_id)
                    .where(TranscriptORM.transcript_id == transcript_row.transcript_id)
                )

                ext += 1

        # Delete object xrefs with no matching xref
        query = (
            select(ObjectXrefCORM.object_xref_id)
            .outerjoin(XrefCORM, XrefCORM.xref_id == ObjectXrefCORM.xref_id)
            .where(XrefCORM.xref_id == None)
        )
        result = core_dbi.execute(query).fetchall()
        object_xref_ids = [row[0] for row in result]

        core_dbi.execute(
            delete(ObjectXrefCORM).where(
                ObjectXrefCORM.object_xref_id.in_(object_xref_ids)
            )
        )

        core_dbi.close()
