import sys
from os import path

from core.input import ServeArgs


def run_server(
    args: ServeArgs,
    base_dir: str,
    nodesdb_f: str,
    pfam_mapping_f: str,
    ipr_mapping_f: str,
    go_mapping_f: str,
    cluster_f: str,
    taxon_idx_mapping_file: str,
    sequence_ids_f: str,
):
    import uvicorn
    from fastapi import FastAPI

    from api.endpoints import router
    from api.sessions import session_manager

    session_manager.base_dir = base_dir
    session_manager.cluster_f = cluster_f
    session_manager.sequence_ids_f = sequence_ids_f
    session_manager.taxon_idx_mapping_file = taxon_idx_mapping_file
    session_manager.nodesdb_f = nodesdb_f
    session_manager.pfam_mapping_f = pfam_mapping_f
    session_manager.ipr_mapping_f = ipr_mapping_f
    session_manager.go_mapping_f = go_mapping_f

    app = FastAPI()

    @app.get("/")
    def hello():
        return {"hi": "hello"}

    app.include_router(router)

    uvicorn.run(app=app, port=args.port)
