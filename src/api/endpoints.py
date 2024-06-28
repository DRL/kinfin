import os
from typing import Dict, List

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.security import APIKeyHeader
from pydantic import BaseModel

from api.sessions import session_manager
from core.input import InputData
from core.results import analyse


class InputSchema(BaseModel):
    config: List[Dict[str, str]]


# X-Session-ID header will be required to access plots/files later
header_scheme = APIKeyHeader(name="x-session-id")

router = APIRouter()


@router.post("/init")
async def initialize(
    input_data: InputSchema,
    background_tasks: BackgroundTasks,
) -> JSONResponse | HTTPException:
    """
    Initialize the analysis process.

    Args:
        input_data (InputSchema): The input data for analysis.
        background_tasks (BackgroundTasks): FastAPI's BackgroundTasks for running analysis asynchronously.

    Returns:
        JSONResponse: A response indicating that the analysis task has been queued.

    Raises:
        HTTPException: If there's an error in the input data or during processing.
    """
    try:
        if not isinstance(input_data.config, list):
            raise HTTPException(
                status_code=400,
                detail="Data must be a list of dictionaries.",
            )

        if not all(isinstance(item, dict) for item in input_data.config):
            raise HTTPException(
                status_code=400,
                detail="Each item in data must be a dictionary.",
            )

        session_id, result_dir = session_manager.new()
        data = InputData(
            nodesdb_f=session_manager.nodesdb_f,
            go_mapping_f=session_manager.go_mapping_f,
            pfam_mapping_f=session_manager.pfam_mapping_f,
            sequence_ids_file=session_manager.sequence_ids_f,
            ipr_mapping_f=session_manager.ipr_mapping_f,
            cluster_file=session_manager.cluster_f,
            config_data=input_data.config,
            taxon_idx_mapping_file=session_manager.taxon_idx_mapping_file,
            output_path=result_dir,
            plot_format="png",  # as we require images
        )

        background_tasks.add_task(analyse, data)

        return JSONResponse(
            content={"detail": "Analysis task has been queued."},
            headers={"X-Session-ID": session_id},
            status_code=202,
        )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")


@router.get("/plot/{plot_type}")
async def get_plot(
    plot_type: str,
    session_id: str = Depends(header_scheme),
) -> FileResponse | HTTPException:
    """
    Retrieve a specific plot type for a given session.

    Args:
        plot_type (str): The type of plot to retrieve.
        session_id (str): The session ID for authentication.

    Returns:
        FileResponse: The requested plot file.

    Raises:
        HTTPException: If the plot type is invalid, session ID is invalid, or the file is not found.
    """
    if plot_type not in ["cluster-size-distribution", "all-rarefaction-curve"]:
        raise HTTPException(status_code=404)

    result_dir = session_manager.get(session_id)
    if not result_dir:
        raise HTTPException(status_code=401, detail="Invalid Session ID provided")

    file_path: str = ""
    match plot_type:
        case "cluster-size-distribution":
            file_path = "cluster_size_distribution.png"
        case "all-rarefaction-curve":
            file_path = "all/all.rarefaction_curve.png"

    file_path = os.path.join(result_dir, file_path)

    if not os.path.exists(file_path):
        raise HTTPException(
            status_code=404,
            detail=f"{plot_type} File Not Found",
        )

    return FileResponse(file_path, media_type="image/png")
