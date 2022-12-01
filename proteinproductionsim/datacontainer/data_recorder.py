"""
================
data_recorder.py
================
"""
import matplotlib.axes
import numpy as np

from ..interface import DataContainer, Environment, Controller
import matplotlib.pyplot as plt
from ..entity.dna_strand import DNAStrand
from ..environment.dna_sim_environment import DNASimEnvironment
from ..variables import data_collection_interval, dt


class DataRecorder(DataContainer):
    """
    This is the basis class for the data recorder. This is a successor class of the DataContainer class.

    :param controller: The controller instance that owns this class. Stored in order to use callback
    """
    def __init__(self, controller: Controller, target):
        super().__init__(controller)
        self.target = target
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
        super().__init__(controller, target)
        self.collection_interval = int(data_collection_interval/dt)
        self._tot_time = total_time + 1
        self._target = target.RNAP_LIST
        self._dt = data_collection_interval
        self._data_list = []
        self._serial_number_list = []
        self._initialized = False
        self._processed_position = None
        self._processed_serial_numbers = None
        self._processed_time_list = None
        self._size = 0

    def log(self, time_index: int):
        # STEP: check if it is time for loading
        if time_index % self.collection_interval == 0:
            # STEP: record the data
            data, serial_number = self._target.get_position_for_recorder()
            self._data_list.append(data)
            self._serial_number_list.append(serial_number)

    def plot(self, axe):
        # STEP: check if the collected data is processed
        if not self._initialized:
            # STEP: we need to convert the raw data into a form which can be plotted
            self._processed_position = []
            self._processed_serial_numbers = []
            self._processed_time_list = []
            cycle_size = len(self._serial_number_list)
            print(cycle_size)
            for i in range(cycle_size):
                true_time = i*self._dt
                for j in range(len(self._data_list[i])):
                    serial_n = self._serial_number_list[i][j]
                    data = self._data_list[i][j]
                    if serial_n not in self._processed_serial_numbers:
                        self._processed_position.append([])
                        self._processed_serial_numbers.append(serial_n)
                        self._processed_time_list.append([])
                        self._size += 1
                    index = self._processed_serial_numbers.index(serial_n)
                    self._processed_position[index].append(data)
                    self._processed_time_list[index].append(true_time)
            self._initialized = True

        axe.set_xlabel('Time [s]')
        axe.set_ylabel('Position [bps]')
        axe.set_title('RNAP Position Plot')
        axe.grid(True)
        for i in range(self._size):
            label = format_n_th_rnap(self._processed_serial_numbers[i])
            axe.plot(self._processed_time_list[i], self._processed_position[i], label=label)
        pass


# Amount Recorder Single-Value Recorder.
class SingleValueRecorder(DataRecorder):
    def __init__(self, controller: Controller, target_function, total_time: int,
                 name_x: str, name_y: str,
                 unit_y: str, unit_x: str,
                 data_format: str = "versus_time"
                 ):
        super().__init__(controller, None)
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
        super().__init__(controller,  target)
        self._target = target.RNAP_LIST
        self._total_time = total_time
        self._name = []
        self._length = 0
        self._loaded = []
        self._detached = []
        self._degrading = []
        self._degraded = []
        self._dt = dt

    def log(self, time_index: int):
        self._loaded.append(self._target.loaded)
        self._detached.append(self._target.detached)
        self._degrading.append(self._target.degrading)
        self._degraded.append(self._target.degraded)

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


class SupercoilingRecorder(DataRecorder):
    """
    This class is used to record supercoiling experienced by the RNAP molecules

    :param controller: the parent Controller Class
    :param target: the target DNAStrand instance
    :param total_time: the integer total time steps
    :param rnap_record_amount: the max amount of rnap to record. default is 5
    """
    def __init__(self, controller, target: DNAStrand, total_time: int, rnap_record_amount: int = 5):
        """
        Constructor method
        """
        super().__init__(controller, target)
        self.collection_interval = int(data_collection_interval / dt)
        self.rnap_record_amount = rnap_record_amount
        self._tot_time = total_time + 1
        self._target = target
        self._dt = data_collection_interval
        self._data_list = []
        self._serial_number_list = []
        self._initialized = False
        self._processed_supercoiling = None
        self._processed_serial_numbers = None
        self._processed_time_list = None
        self._size = 0

    def log(self, time_index: int):
        # STEP: check if it is time for loading
        if time_index % self.collection_interval == 0:
            # STEP: record the data
            data = self.target.phi
            serial_number = self.target.serial_number
            self._data_list.append(data)
            self._serial_number_list.append(serial_number)

    def plot(self, axe):
        # STEP: check if the collected data is processed
        if not self._initialized:
            # STEP: we need to convert the raw data into a form which can be plotted
            self._processed_supercoiling = []
            self._processed_serial_numbers = []
            self._processed_time_list = []
            cycle_size = len(self._serial_number_list)
            print(cycle_size)
            for i in range(cycle_size):
                true_time = i*self._dt
                for j in range(len(self._data_list[i])-1):
                    serial_n = self._serial_number_list[i][j]
                    data = self._data_list[i][j+1]-self._data_list[i][j]
                    if serial_n not in self._processed_serial_numbers:
                        self._processed_supercoiling.append([])
                        self._processed_serial_numbers.append(serial_n)
                        self._processed_time_list.append([])
                        self._size += 1
                    index = self._processed_serial_numbers.index(serial_n)
                    self._processed_supercoiling[index].append(data)
                    self._processed_time_list[index].append(true_time)
            self._initialized = True

        axe.set_xlabel('Time [s]')
        axe.set_ylabel('Supercoiling')
        axe.set_title('RNAP Supercoiling Plot')
        axe.grid(True)
        for i in range(max(self._size-self.rnap_record_amount, 0), self._size):
            label = format_n_th_rnap(self._processed_serial_numbers[i])
            axe.plot(self._processed_time_list[i], self._processed_supercoiling[i], label=label)
        pass


def format_n_th_rnap(number):
    if number % 10 == 0:
        label = f"{number + 1}st RNAP"
    elif number % 10 == 1:
        label = f"{number + 1}nd RNAP"
    elif number % 10 == 2:
        label = f"{number + 1}rd RNAP"
    else:
        label = f"{number + 1}th RNAP"
    return label
