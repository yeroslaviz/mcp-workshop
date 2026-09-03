#!/usr/bin/env python3

import sys
from json import load
from pathlib import Path
from subprocess import run
from typing import Annotated
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from mcp.server import MCPServer
from pydantic import Field


HUMAN_PROTEIN_DATA_DIRECTORY = Path(__file__).parent / "data"
HUMAN_BIOLOGICAL_PROCESS_TERMS = (
    HUMAN_PROTEIN_DATA_DIRECTORY / "9606.protein.enrichment.terms.v12.0.txt"
)
HUMAN_BIOLOGICAL_PROCESS_TERMS_URL = (
    "https://stringdb-downloads.org/download/protein.enrichment.terms.v12.0/"
    "9606.protein.enrichment.terms.v12.0.txt.gz"
)


# Create the MCP server.
mcp = MCPServer("string-workshop")


@mcp.tool(title="STRING: Search human Biological Process terms")
def search_human_biological_process_terms(
    query: Annotated[
        str,
        Field(
            description=(
                "Required. Text to search for in human Gene Ontology Biological Process term "
                "descriptions. Searches case-insensitively. Example: cell cycle"
            )
        ),
    ],
    limit: Annotated[
        int,
        Field(
            description="Optional. Maximum number of matching terms to return. Default: 5.",
            ge=1,
            le=20,
        ),
    ] = 5,
) -> dict:
    """
    Finds human Gene Ontology Biological Process terms whose descriptions match
    the supplied text.

    Each match contains a GO ID and term description. Results are deduplicated
    across STRING proteins.
    """
    print(
        "Tool search_human_biological_process_terms called with parameters: "
        f"query={query}, limit={limit}",
        file=sys.stderr,
    )

    search_text = query.strip()
    if not search_text:
        return {"error": "query must not be empty"}

    if not HUMAN_BIOLOGICAL_PROCESS_TERMS.is_file():
        return {"error": "Human Biological Process term data is unavailable."}

    matches: set[tuple[str, str]] = set()
    more_matches = False

    with open(HUMAN_BIOLOGICAL_PROCESS_TERMS, encoding="utf-8") as term_file:
        for line in term_file:
            _, _, go_id, description = line.rstrip().split("\t", maxsplit=3)
            if search_text.casefold() not in description.casefold():
                continue

            match = (go_id, description)
            if match in matches:
                continue

            if len(matches) >= limit:
                more_matches = True
                break

            matches.add(match)

    return {
        "query": search_text,
        "species": 9606,
        "matches": [
            {"go_id": go_id, "description": description}
            for go_id, description in sorted(matches)
        ],
        "more_matches": more_matches,
    }


@mcp.tool(title="STRING: Get Gene Ontology Biological Process annotations")
def string_biological_process_annotations(
    identifier: Annotated[
        str,
        Field(
            description=(
                "Required. One identifier accepted by STRING. Example: CDK1"
            )
        ),
    ],
    species: Annotated[
        int,
        Field(
            description="Optional. NCBI taxonomy identifier for the proteins. Default: 9606 (human).",
            ge=1,
        ),
    ] = 9606,
    limit: Annotated[
        int,
        Field(
            description=(
                "Optional. Maximum number of Biological Process annotations to return. Default: 10."
            ),
            ge=1,
            le=50,
        ),
    ] = 10,
) -> dict:
    """
    Gets Gene Ontology Biological Process annotations for one protein.
    """
    print(
        "Tool string_biological_process_annotations called with parameters: "
        f"identifier={identifier}, species={species}, limit={limit}",
        file=sys.stderr,
    )

    identifier = identifier.strip()

    if not identifier or "," in identifier or "\n" in identifier:
        return {"error": "Provide exactly one identifier."}

    parameters = urlencode(
        {
            "identifiers": identifier,
            "species": species,
            "caller_identity": "eccb_mcp_workshop",
        }
    )
    url = f"https://string-db.org/api/json/functional_annotation?{parameters}"

    try:
        with urlopen(url, timeout=20) as response:
            annotations = load(response)
    except HTTPError as error:
        return {"error": f"Could not retrieve annotations (HTTP {error.code})."}
    except URLError as error:
        return {"error": "Could not retrieve annotations."}

    biological_process_annotations = [
        annotation
        for annotation in annotations
        if annotation["category"] == "Process"
    ]

    return {
        "input_identifier": identifier,
        "species": species,
        "biological_process_annotations": biological_process_annotations[:limit],
        "more_annotations": len(biological_process_annotations) > limit,
    }


if __name__ == "__main__":
    HUMAN_PROTEIN_DATA_DIRECTORY.mkdir(exist_ok=True)

    if not HUMAN_BIOLOGICAL_PROCESS_TERMS.exists():
        print("Downloading STRING v12.0 human Biological Process terms...", file=sys.stderr)
        command = (
            f"curl -L {HUMAN_BIOLOGICAL_PROCESS_TERMS_URL} | gzip -dc | "
            f"grep 'Biological Process' > {HUMAN_BIOLOGICAL_PROCESS_TERMS}"
        )
        if run(command, shell=True).returncode != 0:
            sys.exit(1)

    if not HUMAN_BIOLOGICAL_PROCESS_TERMS.exists():
        print("Could not download human Biological Process terms.", file=sys.stderr)
        sys.exit(1)

    try:
        mcp.run(
            transport="streamable-http",
            host="0.0.0.0",
            port=8000,
            stateless_http=True,
        )
    except KeyboardInterrupt:
        pass
