from core.input import ServeArgs


def run_server(
    args: ServeArgs,
    results_base_dir: str,
    nodesdb_f: str,
    pfam_mapping_f: str,
    ipr_mapping_f: str,
    go_mapping_f: str,
    cluster_f: str,
    taxon_idx_mapping_file: str,
    sequence_ids_f: str,
) -> None:
    """
    Starts the uvicorn server

    Parameters:
    - args [ServeArgs] : An object containing server configuration arguments, such as the port.
    - results_base_dir [str] : Base directory for storing server results.
    - nodesdb_f [str] : File path to the nodesDB file.
    - pfam_mapping_f [str] : File path to the PFAM mapping file.
    - ipr_mapping_f [str] : File path to the InterPro mapping file.
    - go_mapping_f [str] : File path to the Gene Ontology mapping file.
    - cluster_f [str] : File path to the clustering data file.
    - taxon_idx_mapping_file [str] : File path to the taxon index mapping file.
    - sequence_ids_f [str] : File path to the sequence IDs file.
    """
    import uvicorn
    from fastapi import FastAPI

    from api.endpoints import router
    from api.sessions import session_manager

    session_manager.results_base_dir = results_base_dir
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
