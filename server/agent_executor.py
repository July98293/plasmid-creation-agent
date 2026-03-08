from typing import Optional

from a2a.server.agent_execution import AgentExecutor, RequestContext
from a2a.server.events import EventQueue
from a2a.utils import new_agent_text_message
from openai import AsyncOpenAI


_client: Optional[AsyncOpenAI] = None


def _get_client() -> AsyncOpenAI:
    global _client
    if _client is None:
        _client = AsyncOpenAI()
    return _client


class StablecoinMarketAgent:
    """Stablecoin market simulation agent powered by OpenAI gpt-4.1-nano.

    Simulates short-term market conditions for major USD stablecoins
    (USDC, USDT, DAI), including prices, yields, and depeg risks,
    in response to a portfolio manager's actions or questions.
    """

    def __init__(self) -> None:
        self._client = _get_client()

    async def invoke(self, user_input: str) -> str:
        prompt = (
            "You are a STABLECOIN MARKET SIMULATOR for USDC, USDT, and DAI.\n"
            "Time horizon: short-term (next few hours to 1 day).\n"
            "The portfolio manager will describe its current holdings, planned trades,\n"
            "or questions about the market. Your job is to simulate a realistic but\n"
            "fictional market response for that horizon.\n\n"
            "Always include:\n"
            "- Approximate prices for USDC, USDT, and DAI vs USD (e.g. 1.00, 0.998, 1.002).\n"
            "- Any notable yield opportunities (e.g. lending/LP yields) and where.\n"
            "- Any depeg or regulatory risk events that might affect the portfolio.\n"
            "- A brief explanation of why the market moved that way.\n\n"
            f"Portfolio manager message:\n{user_input or 'Describe a baseline stable, low-volatility market.'}"
        )

        response = await self._client.responses.create(
            model="gpt-4.1-nano",
            input=prompt,
        )

        text = getattr(response, "output_text", None)
        if text is None:
            try:
                first_item = response.output[0]
                content = first_item.content[0]
                # Newer SDKs wrap text in .text.value
                if getattr(content, "text", None) is not None:
                    text = getattr(content.text, "value", None) or str(content.text)
                else:
                    text = str(content)
            except Exception:  # pragma: no cover - defensive fallback
                text = str(response)

        return text


class StablecoinMarketAgentExecutor(AgentExecutor):
    """AgentExecutor that wraps the StablecoinMarketAgent."""

    def __init__(self) -> None:
        self._agent = StablecoinMarketAgent()

    async def execute(
        self,
        context: RequestContext,
        event_queue: EventQueue,
    ) -> None:
        user_input = context.get_user_input()
        result_text = await self._agent.invoke(user_input)
        await event_queue.enqueue_event(new_agent_text_message(result_text))

    async def cancel(
        self,
        context: RequestContext,
        event_queue: EventQueue,
    ) -> None:
        # No special cancel behaviour for this simple agent.
        return None

