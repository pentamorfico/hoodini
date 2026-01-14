from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Literal

import polars as pl

Mode = Literal["ipg", "nuc2ass", "nucsum", "asssum"]


class IPGXMLFile:
    def __init__(self, file_path: str | Path):
        """Parse an XML file from NCBI IPG/eutils responses."""

        self.file_path = Path(file_path)
        self.xml_tree = ET.parse(self.file_path)
        self.root = self.xml_tree.getroot()

    def to_dict(self, mode: Mode) -> dict[str, Any]:
        """Return the parsed XML content as a nested dictionary."""

        if mode not in ("ipg", "nuc2ass", "nucsum", "asssum"):
            raise ValueError(f"Unsupported mode '{mode}'")

        parsed_dict: dict[str, Any] = {}

        if mode == "ipg":
            for ipg_report in self.root.findall("IPGReport"):
                ipg_id = ipg_report.get("ipg")
                product_acc = ipg_report.get("product_acc")

                product = ipg_report.find("Product")
                product_details = product.attrib if product is not None else {}

                proteins = []
                for protein in ipg_report.findall(".//Protein"):
                    protein_accver = protein.get("accver")
                    protein_details = {k: v for k, v in protein.attrib.items() if k != "accver"}

                    cds_list = []
                    for cds in protein.findall(".//CDS"):
                        cds_list.append(cds.attrib)

                    proteins.append(
                        {
                            "protein_accver": protein_accver,
                            "protein_details": protein_details,
                            "cds_list": cds_list,
                        }
                    )

                parsed_dict[ipg_id] = {
                    "product_acc": product_acc,
                    "product_details": product_details,
                    "proteins": proteins,
                }

        elif mode == "nuc2ass":
            for link_set in self.root.findall("LinkSet"):
                db_from = link_set.findtext("DbFrom")
                id_list = [id_tag.text for id_tag in link_set.findall(".//IdList/Id")]

                link_set_db = link_set.find("LinkSetDb")
                if link_set_db is not None:
                    db_to = link_set_db.findtext("DbTo")
                    link_name = link_set_db.findtext("LinkName")
                    linked_ids = [link.findtext("Id") for link in link_set_db.findall(".//Link")]
                else:
                    db_to = link_name = None
                    linked_ids = []

                parsed_dict[tuple(id_list)] = {
                    "db_from": db_from,
                    "db_to": db_to,
                    "link_name": link_name,
                    "linked_ids": linked_ids,
                }

        elif mode == "nucsum":
            for doc_sum in self.root.findall("DocSum"):
                doc_id = doc_sum.findtext("Id")
                doc_items = {item.get("Name"): item.text for item in doc_sum.findall("Item")}
                parsed_dict[doc_id] = doc_items

        elif mode == "asssum":
            for doc_summary in self.root.findall(".//DocumentSummary"):
                uid = doc_summary.get("uid")
                entry = {child.tag: child.text for child in doc_summary}
                parsed_dict[uid] = entry

        return parsed_dict

    def to_dataframe(self, mode: Mode) -> pl.DataFrame:
        """Flatten the parsed XML content into a Polars DataFrame."""

        parsed_dict = self.to_dict(mode)
        flattened_data: list[dict[str, Any]] = []

        if mode == "ipg":
            for ipg_id, details in parsed_dict.items():
                base_info = {"ipg_id": ipg_id, "product_acc": details["product_acc"]}
                base_info.update(details["product_details"])

                for protein in details["proteins"]:
                    protein_info = protein["protein_details"].copy()
                    protein_info["protein_accver"] = protein["protein_accver"]

                    if protein["cds_list"]:
                        for cds in protein["cds_list"]:
                            flattened_data.append({**base_info, **protein_info, **cds})
                    else:
                        flattened_data.append({**base_info, **protein_info})

        elif mode == "nuc2ass":
            for ids, details in parsed_dict.items():
                id_list_str = ",".join(ids)
                for linked_id in details["linked_ids"]:
                    flattened_data.append(
                        {
                            "db_from": details["db_from"],
                            "id_list": id_list_str,
                            "db_to": details["db_to"],
                            "link_name": details["link_name"],
                            "linked_id": linked_id,
                        }
                    )

        elif mode == "nucsum":
            for doc_id, doc_items in parsed_dict.items():
                row = {"doc_id": doc_id}
                row.update(doc_items)
                flattened_data.append(row)

        elif mode == "asssum":
            for uid, details in parsed_dict.items():
                row = {"uid": uid}
                row.update(details)
                flattened_data.append(row)

        return pl.DataFrame(flattened_data)
