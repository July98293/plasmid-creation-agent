from __future__ import annotations

import json
from typing import Any, Dict, Optional

import httpx

from logging_utils import get_conversation_logger


async def fetch_json(
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    label: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    HTTP GET wrapper that logs request/response pairs to the shared
    conversation log so that all "tool calls" are inspectable.
    """
    conv_logger = get_conversation_logger()
    tag = label or "HTTP_JSON"

    conv_logger.info(
        "[%s] REQUEST url=%s params=%s headers=%s",
        tag,
        url,
        json.dumps(params or {}),
        json.dumps(headers or {}),
    )

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            resp = await client.get(url, params=params, headers=headers)
    except Exception as exc:  # noqa: BLE001
        conv_logger.info(
            "[%s] ERROR url=%s error=%s",
            tag,
            url,
            repr(exc),
        )
        return None

    body_text: str
    try:
        data = resp.json()
        body_text = json.dumps(data)[:2000]
    except Exception:  # noqa: BLE001
        data = None
        body_text = (resp.text or "")[:2000]

    conv_logger.info(
        "[%s] RESPONSE url=%s status=%s body=%s",
        tag,
        url,
        resp.status_code,
        body_text,
    )

    if not resp.is_success:
        return None

    return data


async def fetch_text(
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    headers: Optional[Dict[str, str]] = None,
    label: Optional[str] = None,
) -> Optional[str]:
    """
    HTTP GET wrapper for text endpoints (e.g. FASTA) with the same logging
    guarantees as fetch_json.
    """
    conv_logger = get_conversation_logger()
    tag = label or "HTTP_TEXT"

    conv_logger.info(
        "[%s] REQUEST url=%s params=%s headers=%s",
        tag,
        url,
        json.dumps(params or {}),
        json.dumps(headers or {}),
    )

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            resp = await client.get(url, params=params, headers=headers)
    except Exception as exc:  # noqa: BLE001
        conv_logger.info(
            "[%s] ERROR url=%s error=%s",
            tag,
            url,
            repr(exc),
        )
        return None

    text = (resp.text or "")[:2000]
    conv_logger.info(
        "[%s] RESPONSE url=%s status=%s body=%s",
        tag,
        url,
        resp.status_code,
        text,
    )

    if not resp.is_success:
        return None

    return resp.text

