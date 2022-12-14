"""
=====================
dna_sim_controller.py
=====================

This file defines the main controller for the simulation.
"""

from proteinproductionsim.interface import Controller, DataContainer
from ..environment.dna_sim_environment import DNASimEnvironment
from ..variables import stage_per_collection, dt, total_time, scaling, t_on
from ..datacontainer.data_recorder import RNAPPositionRecorder, SingleValueRecorder, FiveThreeRecorder, \
    SupercoilingRecorder


class RecordConfig:
    """
    This class is used to pass and store data recording setting to the simulator.
    """
    def __init__(self, controller=None, record_rnap_position: bool = False, record_rnap_amount: bool = False,
                 record_rnap_state: bool = True, record_processing_time: bool = True,
                 record_protein_amount: bool = True, record_protein_production: bool = True,
                 record_finish_time: bool = True, record_five_three: bool = False, record_supercoiling: bool = False,
                 log_data: bool = True, show_progress_bar: bool = True):
        self.parent = controller
        self.record_rnap_position = record_rnap_position
        self.record_rnap_amount = record_rnap_amount
        self.record_rnap_state = record_rnap_state
        self.record_processing_time = record_processing_time
        self.record_protein_amount = record_protein_amount
        self.record_protein_production = record_protein_production
        self.record_finish_time = record_finish_time
        self.record_five_three = record_five_three
        self.log_data = log_data
        self.show_progress_bar = show_progress_bar
        self.record_supercoiling = record_supercoiling


class RunConfig:
    """
    This class is used to pass
    """
    def __init__(self,
                 rnap_loading_rate: float = 0.127):
        self.rnap_loading_rate = rnap_loading_rate
        self.config = {}

    def get(self, key):
        return self.config[key]


class DNASimController(Controller):
    """
    This is the central controller for the DNA simulation.
    """
    def __init__(self, rnap_loading_rate: float, record_config: RecordConfig = RecordConfig(), **kwargs):
        super().__init__()
        # Setup
        self.time_index = 0
        self.env = DNASimEnvironment(controller=self, rnap_loading_rate=rnap_loading_rate,
                                     if_storing_supercoiling_value=record_config.record_supercoiling,
                                     **kwargs)
        self.total_time = scaling(total_time)
        self.dt = dt
        self.stage_per_collection = stage_per_collection

        # Data Recording Setup
        self.record_config = record_config

        # Initialize Data Recorders
        self.data_recorder = {}
        pass

    def start(self):
        self.init()
        for i in range(self.total_time):
            self.env.step(time_index=self.time_index)
            self._log()
            self.time_index += 1

        pass

    def init(self):
        # adding all the data recorder according to the config
        if self.record_config.record_rnap_position:
            self.data_recorder["position"] = RNAPPositionRecorder(self, self.env.dna, self.total_time)
        if self.record_config.record_protein_amount:
            self.data_recorder["protein amount"] = SingleValueRecorder(self, self.get_protein_amount, self.total_time,
                                                                       name_x="Time", name_y="Protein Amount",
                                                                       unit_x="s", unit_y="")
        if self.record_config.record_five_three:
            self.data_recorder["five and three"] = FiveThreeRecorder(self, self.env.dna, self.total_time)
        if self.record_config.record_supercoiling:
            self.data_recorder["supercoiling"] = SupercoilingRecorder(self, self.env.dna, self.total_time)

        self.env.init()
        pass

    def call_back(self, option, data):
        pass

    def _log(self):
        for key in self.data_recorder:
            self.data_recorder[key].log(self.time_index)
        pass

    def get_data(self, name):
        if name in self.data_recorder:
            return self.data_recorder[name]
        return None

    def get_protein_amount(self):
        return self.env.total_prot

    def get_five_three(self):
        if self.record_config.record_five_three:
            return self.data_recorder["five and three"].get_five_six()
        else:
            return None
