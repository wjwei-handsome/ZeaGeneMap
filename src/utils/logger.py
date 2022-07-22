import logging
from rich.logging import RichHandler
from rich import print
from rich.console import Console

console = Console()
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[RichHandler(rich_tracebacks=True)])
logger = logging.getLogger("GeneMap")
pretty_print = print.__call__
