from typing import Dict, List

from fastapi import APIRouter
from pydantic import BaseModel

from api.sessions import session_manager
from core.input import InputData
from core.results import analyse


class InputSchema(BaseModel):
    data: List[Dict[str, str]]


router = APIRouter()


@router.get("/init")
async def initialize(input_data: InputSchema):
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
    )
    analyse(data)
    return {"message": session_id}
