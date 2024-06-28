import os
import shutil
import signal
import threading
import time
from datetime import datetime, timedelta
from typing import Dict, Optional, Tuple
from uuid import uuid4

from core.utils import logger


class Session:
    # TODO : Update to use Query Session instead of User Sessions
    """
    Represents a user session.

    Attributes:
        session_id (str): Unique identifier for the session.
        result_path (str): Path where session results are stored.
        last_activity (datetime): Timestamp of the last activity in the session.
    """

    def __init__(self, results_base_dir: str) -> None:
        """
        Initialize a new session.

        Args:
            results_base_dir (str): Base directory where session results will be stored.
        """
        self.session_id = uuid4().hex

        self.result_path = os.path.join(results_base_dir, self.session_id)
        os.makedirs(self.result_path, exist_ok=True)
        self.last_activity = datetime.now()

    def update_activity(self) -> None:
        """Update the timestamp of the last activity in the session"""
        self.last_activity = datetime.now()

    def is_expired(self) -> bool:
        """
        Check if the session has expired based on a threshold.

        Returns:
            bool: True if the session is expired, False otherwise.
        """
        threshold = os.getenv("SESSION_INACTIVITY_THRESHOLD")
        default_threshold_hours = 24
        try:
            if threshold is None:
                logger.warning(
                    "[WARN] - SESSION_INACTIVITY_THRESHOLD environment variable not set. Defaulting to 24 hours."
                )
                threshold_hours = default_threshold_hours
            else:
                threshold_hours = float(threshold)
        except ValueError:
            logger.error(
                "[ERROR] - Invalid value for SESSION_INACTIVITY_THRESHOLD environment variable. Defaulting to 24 hours."
            )
            threshold_hours = default_threshold_hours

        return datetime.now() > self.last_activity + timedelta(hours=threshold_hours)


class SessionManager:
    """
    Manages user sessions and their lifecycle.

    Attributes:
        results_base_dir (str): Base directory where session results are stored.
        sessions (Dict[str, Session]): Dictionary to store active sessions.
        cleanup_thread (Thread): Thread for periodic cleanup of expired sessions.
    """

    def __init__(self) -> None:
        """Initialize the SessionManager"""
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
        """
        Create a new session and return its ID and result path.

        Returns:
            Tuple[str, str]: ID and result path of the new session.
        """
        session = Session(self.results_base_dir)
        self.sessions[session.session_id] = session
        return session.session_id, session.result_path

    def get(self, session_id) -> Optional[str]:
        """
        Retrieve the result path for a given session ID.

        Args:
            session_id (str): ID of the session to retrieve.

        Returns:
            Optional[str]: Result path if session is active and not expired, None otherwise.
        """
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
        """
        Remove a session by its ID.

        Args:
            session_id (str): ID of the session to remove.

        Returns:
            bool: True if session was successfully removed, False otherwise.
        """
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
        """Remove all expired sessions"""
        expired_sessions = [
            session_id
            for session_id, session in self.sessions.items()
            if session.is_expired()
        ]
        for session_id in expired_sessions:
            self.remove(session_id)

    def clear_all_sessions(self) -> None:
        """Remove all sessions and associated result directories"""
        shutil.rmtree(self.results_base_dir)
        self.sessions = {}

    def cleanup_loop(self) -> None:
        """Periodically clean up expired sessions"""
        while True:
            self.clear_expired_sessions()
            time.sleep(10)

    def __exit__(self, signum, frame) -> None:
        """Cleanup all sessions when exiting due to signal"""
        self.clear_all_sessions()
        exit(0)


session_manager = SessionManager()

# Delete stored results on server shutdown
signal.signal(signal.SIGINT, session_manager.__exit__)
signal.signal(signal.SIGTERM, session_manager.__exit__)
