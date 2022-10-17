from interface import Controller

from interface import Controller
from controller.dna_sim_controller import DNASimController


class MultiSampleController(Controller):
    """
    This is the controller for large sample simulation. This will directly control other single-sample controller.
    """
    def __init__(self, run_function, ):
        super().__init__()

        pass

    def init(self):
        pass

    def start(self):
        pass

    def call_back(self, option, data):
        pass