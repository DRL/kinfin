import os
import shutil
import signal
import threading
import time
from datetime import datetime, timedelta
from typing import Dict, Tuple
from uuid import uuid4


class Session:
    def __init__(self, base_dir: str) -> None:
        self.session_id = uuid4().hex

        self.session_path = os.path.join(base_dir, f"sessions/{self.session_id}")
        os.makedirs(self.session_path, exist_ok=True)
        self.last_activity = datetime.now()

    def update_activity(self) -> None:
        self.last_activity = datetime.now()

    def is_expired(self) -> bool:
        return datetime.now() > self.last_activity + timedelta(hours=3)


class SessionManager:
    def __init__(self) -> None:
        self.base_dir = ""
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
        session = Session(self.base_dir)
        self.sessions[session.session_id] = session
        return session.session_id, session.session_path

    def get(self, id) -> Session | None:
        session = self.sessions.get(id)

        if not session:
            return None

        expired = session.is_expired()

        if not expired:
            session.update_activity()
            return session
        else:
            self.remove(id)
            return None

    def remove(self, id) -> bool:
        shutil.rmtree(self.sessions[id].session_path)
        session = self.sessions.pop(id, None)
        if session:
            try:
                os.rmdir(f"sessions/{id}")
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
        for session_id in list(self.sessions.keys()):
            self.remove(session_id)

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
