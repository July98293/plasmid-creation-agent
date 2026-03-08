import argparse
import logging
import os

from dotenv import load_dotenv

from a2a.server.apps import A2AStarletteApplication
from a2a.server.request_handlers import DefaultRequestHandler
from a2a.server.tasks import InMemoryTaskStore
from a2a.types import AgentCapabilities, AgentCard, AgentSkill

from .agent_executor import StablecoinMarketAgentExecutor


def build_agent_card(port: int) -> AgentCard:
    skill = AgentSkill(
        id="simulate_stablecoin_market",
        name="Simulate stablecoin market",
        description=(
            "Simulates short-term market conditions for USDC, USDT, and DAI, "
            "including prices, yields, and depeg risks."
        ),
        tags=["stablecoin", "market", "simulation"],
        examples=[
            "Given my portfolio in USDC, USDT, and DAI, simulate the next 24 hours.",
            "Simulate a mild depeg event for USDT and its impact on my portfolio.",
        ],
    )

    return AgentCard(
        name="Stablecoin Market Agent",
        description=(
            "Stablecoin market simulator for USDC, USDT, and DAI over the next few hours to 1 day."
        ),
        url=f"http://127.0.0.1:{port}/",
        version="1.0.0",
        default_input_modes=["text"],
        default_output_modes=["text"],
        capabilities=AgentCapabilities(streaming=True),
        skills=[skill],
    )


def main() -> None:
    load_dotenv()
    load_dotenv("../.env")

    parser = argparse.ArgumentParser(
        description="Run an A2A JSON-RPC server using the a2a-python SDK.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=9090,
        help="Port for the A2A server to listen on.",
    )
    args = parser.parse_args()

    # Mirror the Go server's environment variable mapping so existing
    # OPENAI_API_KEY configs continue to work if you later plug in an LLM.
    if not os.getenv("OPENAIAPIKEY") and os.getenv("OPENAI_API_KEY"):
        os.environ["OPENAIAPIKEY"] = os.environ["OPENAI_API_KEY"]

    logging.basicConfig(level=logging.INFO)

    if not os.getenv("OPENAIAPIKEY"):
        logging.warning(
            "[server] OPENAIAPIKEY is not set; stablecoin market agent will not call external LLMs.",
        )
    else:
        logging.info("[server] OPENAIAPIKEY is set (value hidden)")

    agent_card = build_agent_card(args.port)

    executor = StablecoinMarketAgentExecutor()
    request_handler = DefaultRequestHandler(
        agent_executor=executor,
        task_store=InMemoryTaskStore(),
    )

    app = A2AStarletteApplication(
        agent_card=agent_card,
        http_handler=request_handler,
        extended_agent_card=agent_card,
    )

    import uvicorn

    logging.info(
        "Starting A2A JSON-RPC server on http://127.0.0.1:%d",
        args.port,
    )

    uvicorn.run(app.build(), host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()

