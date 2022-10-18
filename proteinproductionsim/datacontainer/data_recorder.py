"""
================
data_recorder.py
================
"""
import numpy as np

from ..interface import DataContainer, Environment, Controller
import matplotlib.pyplot as plt
from ..entity.dna_strand import DNAStrand
from ..environment.dna_sim_environment import DNASimEnvironment
from ..variables import data_collection_interval, dt


class DataRecorder(DataContainer):
    """
    This is the basis class for the data recorder.
    """
    def __init__(self, controller: Controller, collection_interval):
        super().__init__(controller)
        self.collection_interval = collection_interval
        pass

    def init(self):
        pass

    def log(self, time_index):
        pass

    def get(self):
        pass

    def plot(self, fig):
        pass

    def store(self, path):
        pass


class RNAPPositionRecorder(DataRecorder):
    def __init__(self, controller, target: DNAStrand, total_time: int):
        super().__init__(controller, int(data_collection_interval/dt))
        self._tot_time = total_time + 1
        self._dna = target
        self._dt = dt
        self._data = []
        self._length = 0
        self._rnap_loaded = 0
        self._rnap_detached = 0
        self._start_end = []

    def log(self, time_index: int):
        # STEP: check RNAP amount through loaded and attached
        if time_index % self.collection_interval == 0:
            # STEP: check if detached increase
            if self._dna.detached > self._rnap_detached:
                self._start_end[self._rnap_detached][1] = time_index - 1
                self._rnap_detached += 1
            # STEP: check if loaded increase
            if self._dna.loaded > self._rnap_loaded:
                self._start_end.append([time_index, -1])
                self._rnap_loaded += 1
                self._data.append(np.zeros(self._tot_time))

        # STEP: now we track the position.
        for i in range(self._rnap_loaded):
            if i < self._rnap_detached:
                continue
            self._data[i][self._length] = self._dna.RNAP_LIST[i].position

        self._length += 1
        pass

    def plot(self, axe):
        axe.set_xlabel('Time [s]')
        axe.set_ylabel('Position [bps]')
        axe.set_title('RNAP Position Plot')
        axe.grid(True)
        for i in range(self._rnap_loaded):
            start_index = self._start_end[i][0]
            stop_index = self._start_end[i][1]
            if stop_index < 0:
                stop_index = self._tot_time
            time_list = np.arange(start=start_index, stop=stop_index,
                                  step=1, dtype=float)
            axe.plot(time_list*self._dt, self._data[i][start_index:stop_index])
        #print([self._rnap_loaded, self._rnap_detached])
        pass


# Amount Recorder Single-Value Recorder.
class SingleValueRecorder(DataRecorder):
    def __init__(self, controller: Controller, target_function, total_time: int,
                 name_x: str, name_y: str,
                 unit_y: str, unit_x: str,
                 data_format: str = "versus_time"
                 ):
        super().__init__(controller, int(data_collection_interval/dt))
        self._tot_time = total_time
        self._targe_function = target_function
        self._data_format = data_format
        self._name = [name_x, name_y]
        self._unit = [unit_x, unit_y]
        self._data = np.zeros(self._tot_time)

    def log(self, time_index: int) -> None:
        self._data[time_index] = self._targe_function()

    def get(self):
        return self._data

    def plot(self, axe):
        x_label = f"{self._name[0]} [{self._unit[0]}]"
        if self._unit[1] == "":
            y_label = f"{self._name[1]}"
        else:
            y_label = f"{self._name[1]} [{self._unit[1]}]"
        axe.set_xlabel(x_label)
        axe.set_ylabel(y_label)
        axe.set_title(f'{self._name[0]} versus {self._name[1]} Plot')
        axe.grid(True)
        time_list = np.arange(start=0, stop=self._tot_time, step=1, dtype=float)
        axe.plot(time_list*dt, self._data)
        pass


# Multi-Value Recorder
class FiveThreeRecorder(DataRecorder):
    def __init__(self, controller, target: DNAStrand, total_time):
        super().__init__(controller,  int(data_collection_interval/dt))
        self._dna = target
        self._total_time = total_time
        self._name = []
        self._length = 0
        self._loaded = []
        self._detached = []
        self._degrading = []
        self._degraded = []
        self._dt = dt

    def log(self, time_index: int):
        self._loaded.append(self._dna.loaded)
        self._detached.append(self._dna.detached)
        self._degrading.append(self._dna.degrading)
        self._degraded.append(self._dna.degraded)

        self._length += 1
        return 0

    def get_five_six(self):
        five = np.array(self._loaded) - np.array(self._degrading)
        three = np.array(self._detached) - np.array(self._degraded)
        return [five, three]

    def plot(self, axe):
        five, three = self.get_five_six()
        axe.set_xlabel('Time [s]')
        axe.set_ylabel('Amount')
        axe.set_title('5 and 3 mRNA amount versus time Plot')
        axe.grid(True)
        time_list = np.arange(start=0, stop=self._total_time,
                              step=1, dtype=float) * self._dt
        axe.plot(time_list, five, label='Five End')
        axe.plot(time_list, three, label='Three End')
        # print([self._rnap_loaded, self._rnap_detached])
        pass
