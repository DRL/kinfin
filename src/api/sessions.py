import os
import shutil
import signal
import threading
import time
from datetime import datetime, timedelta
from typing import Dict, Optional, Tuple
from uuid import uuid4


class Session:
    def __init__(self, results_base_dir: str) -> None:
        self.session_id = uuid4().hex

        self.result_path = os.path.join(results_base_dir, self.session_id)
        os.makedirs(self.result_path, exist_ok=True)
        self.last_activity = datetime.now()

    def update_activity(self) -> None:
        self.last_activity = datetime.now()

    def is_expired(self) -> bool:
        return datetime.now() > self.last_activity + timedelta(hours=3)


class SessionManager:
    def __init__(self) -> None:
        self.results_base_dir = ""
        self.cluster_f = ""
        self.sequence_ids_f = ""
        self.taxon_idx_mapping_file = ""
        self.nodesdb_f = ""
        self.pfam_mapping_f = ""
        self.ipr_mapping_f = ""
        self.go_mapping_f = ""

        self.sessions: Dict[str, Session] = {}
        self.cleanup_thread = threading.Thread(target=self.cleanup_loop, daemon=True)
        self.cleanup_thread.start()

    def new(self) -> Tuple[str, str]:
        session = Session(self.results_base_dir)
        self.sessions[session.session_id] = session
        return session.session_id, session.result_path

    def get(self, session_id) -> Optional[str]:
        session = self.sessions.get(session_id)

        if not session:
            return None

        expired = session.is_expired()

        if not expired:
            session.update_activity()
            return session.result_path
        else:
            self.remove(session_id)
            return None

    def remove(self, session_id) -> bool:
        shutil.rmtree(self.sessions[session_id].result_path)
        session = self.sessions.pop(session_id, None)
        if session:
            try:
                os.rmdir(session.result_path)
            except FileNotFoundError:
                pass
            return True
        else:
            return False

    def clear_expired_sessions(self) -> None:
        expired_sessions = [
            session_id
            for session_id, session in self.sessions.items()
            if session.is_expired()
        ]
        for session_id in expired_sessions:
            self.remove(session_id)

    def clear_all_sessions(self):
        shutil.rmtree(self.results_base_dir)

    def cleanup_loop(self):
        while True:
            self.clear_expired_sessions()
            time.sleep(10)

    def __exit__(self, signum, frame):
        self.clear_all_sessions()
        exit(0)


session_manager = SessionManager()
signal.signal(signal.SIGINT, session_manager.__exit__)
signal.signal(signal.SIGTERM, session_manager.__exit__)
