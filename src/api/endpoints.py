import os
from typing import Dict, List

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel

from api.sessions import session_manager
from core.input import InputData
from core.results import analyse


class InputSchema(BaseModel):
    data: List[Dict[str, str]]


router = APIRouter()


@router.post("/init")
async def initialize(input_data: InputSchema):
    try:
        if not isinstance(input_data.data, list):
            raise HTTPException(
                status_code=400,
                detail="Data must be a list of dictionaries.",
            )
        for item in input_data.data:
            if not isinstance(item, dict):
                raise HTTPException(
                    status_code=400,
                    detail="Each item in data must be a dictionary.",
                )
        session_id, session_path = session_manager.new()
        data = InputData(
            nodesdb_f=session_manager.nodesdb_f,
            go_mapping_f=session_manager.go_mapping_f,
            pfam_mapping_f=session_manager.pfam_mapping_f,
            sequence_ids_file=session_manager.sequence_ids_f,
            ipr_mapping_f=session_manager.ipr_mapping_f,
            cluster_file=session_manager.cluster_f,
            config_data=input_data.data,
            taxon_idx_mapping_file=session_manager.taxon_idx_mapping_file,
            output_path=session_path,
            plot_format="png",  # as we require images
        )
        analyse(data)
        file_path = os.path.join(session_path, "cluster_size_distribution.png")
        if not os.path.exists(file_path):
            raise HTTPException(
                status_code=500,
                detail="cluster_size_distribution.png could not be generated",
            )

        return FileResponse(
            file_path,
            media_type="image/png",
            headers={"X-Session-ID": session_id},
        )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")
