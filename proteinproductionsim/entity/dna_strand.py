"""
============
dna_strand.py
============

This file contains the DNAStrand class.
"""
from numpy import ndarray

from proteinproductionsim.interface import Entity
from proteinproductionsim.variables import length, scaling, dt, total_time, v_0, pauseDuration, pauseSite, RNAP_size, \
    t_on, tau_0, tau_c, stalling_supercoiling

from proteinproductionsim.helper.loading_list import LoadingList
from proteinproductionsim.helper.supercoilling import s, n_dependence_cubic_3
from proteinproductionsim.helper.general import if_out_of_interval
from proteinproductionsim.entity.rnap import RNAP
import numpy as np


class RNAPList:
    def __init__(self, dna):
        # basic initiation
        self.dna = dna  # keep the parent instance
        self.attached_rnap_list: list[RNAP] = []  # this list stores the attached RNAPs instances
        self.detached_rnap_list: list[RNAP] = []  # this list stores the detached RNAPs instances
        self.inert_rnap_list: list[RNAP] = []  # this list stores all the deactivated RNAPs instances

        # variables related to the status of its children mRNAs.
        self.loaded = 0  # number of RNAPs that have been loaded.
        self.detached = 0  # number of RNAPs that have been detached.
        self.attached = 0  # number of RNAPs that are attached.
        self.degrading = 0  # number of RNAPs that are degrading.
        self.degraded = 0  # number of RNAPs that are already degraded.
        self.interrupted = 0  # number of RNAPs that interrupted.

        # adaptive supercoiling
        self.r_ref: list[float] = []  # reference position of active RNAPs for adaptive supercoiling
        self.flag_r_ref: list[bool] = []  # boolean for the whether the r_ref is used

    def init(self):
        return

    def transfer_element(self, old_list, new_list, element):
        if old_list is self.attached_rnap_list:
            index = self.attached_rnap_list.index(element)
            self.r_ref.pop(index)
            self.flag_r_ref.pop(index)
        old_list.remove(element)
        new_list.append(element)
        return True

    def if_loading_site_clean(self):
        if self.loaded == 0:
            return True
        if_clean = True
        last_attached_rnap: RNAP = self.get_rear_attached_rnap()
        if (last_attached_rnap.position - RNAP_size) < 0:
            if_clean = False
        return if_clean

    def attach_rnap(self, **kwargs):
        serial_n = self.loaded
        self.attached_rnap_list.append(RNAP(self, serial_n=serial_n, **kwargs))
        self.loaded += 1
        self.attached += 1
        self.r_ref.append(0)
        self.flag_r_ref.append(False)
        pass

    def get_attached_serial_number(self) -> list[int]:
        attached_rnap_serial_number: list[int] = []
        for i in range(len(self.attached_rnap_list)):
            attached_rnap_serial_number.append(self.attached_rnap_list[i].serial_number)
        return attached_rnap_serial_number

    def get_detached_serial_number(self) -> list[int]:
        detached_rnap_serial_number: list[int] = []
        for i in range(len(self.detached_rnap_list)):
            detached_rnap_serial_number.append(self.detached_rnap_list[i].serial_number)
        return detached_rnap_serial_number

    def get_position_for_all_attached_rnap(self) -> tuple[list[float], list[int]]:
        attached_rnap_position: list[float] = []
        attached_rnap_serial_number: list[int] = []
        for i in range(len(self.attached_rnap_list)):
            rnap = self.attached_rnap_list[i]
            attached_rnap_serial_number.append(rnap.serial_number)
            attached_rnap_position.append(rnap.position)

        return attached_rnap_position, attached_rnap_serial_number

    def get_supercoiling_ref(self):
        return self.r_ref, self.flag_r_ref

    def get_attached_rnap(self, serial_number):
        for rnap in self.attached_rnap_list:
            if rnap.serial_number == serial_number:
                return rnap

    def get_detached_rnap(self, serial_number):
        for rnap in self.detached_rnap_list:
            if rnap.serial_number == serial_number:
                return rnap

    def step(self, time_index, stepping: list[float], serial_number_list: list[int]) -> int:
        prot = 0
        for i in range(len(serial_number_list)):
            serial_number = serial_number_list[i]
            rnap = self.get_attached_rnap(serial_number)
            prot += rnap.step(time_index, stepping[i])
        detached_rnap_serial_number_list = self.get_detached_serial_number()
        for i in range(len(detached_rnap_serial_number_list)):
            serial_number = detached_rnap_serial_number_list[i]
            rnap = self.get_detached_rnap(serial_number)
            prot += rnap.step(time_index, 0.0)
        return prot

    def call_back(self, operation: str, entity: RNAP):
        match operation:
            case "attached":
                pass
            case "detached":
                self.attached -= 1
                self.detached += 1
                self.transfer_element(self.attached_rnap_list, self.detached_rnap_list, entity)
            case "degrading":
                self.degrading += 1
            case "degraded":
                self.degraded += 1
                self.transfer_element(self.detached_rnap_list, self.inert_rnap_list, entity)
            case "interrupted":
                self.attached -= 1
                self.interrupted += 1
                self.process_rnap_interruption(entity)
            case "high_supercoiling_fall_off":
                self.attached -= 1
                self.interrupted += 1
                self.process_rnap_interruption(entity)

    def process_rnap_detachment(self, entity: RNAP):
        idx = self.attached_rnap_list.index(entity)
        self.transfer_element(self.attached_rnap_list, self.detached_rnap_list, entity)
        self.r_ref.pop(idx)
        self.flag_r_ref.pop(idx)

    def process_rnap_interruption(self, entity: RNAP):

        if entity in self.attached_rnap_list:
            self.accumulate_r_ref_to_the_front_rnap(entity)
            self.transfer_element(self.attached_rnap_list, self.inert_rnap_list, entity)
        elif entity in self.detached_rnap_list:
            self.transfer_element(self.detached_rnap_list, self.detached_rnap_list, entity)

    def get_rear_attached_rnap(self) -> RNAP:
        return self.attached_rnap_list[-1]

    def get_position_for_recorder(self):
        data = np.zeros(self.attached)
        serial_number = np.zeros(self.attached, dtype=int)
        for i in range(len(self.attached_rnap_list)):
            data[i] = self.attached_rnap_list[i].position
            serial_number[i] = self.attached_rnap_list[i].serial_number
        return data, serial_number

    def accumulate_r_ref_to_the_front_rnap(self, entity: RNAP):
        # STEP: get RNAP index
        idx = self.attached_rnap_list.index(entity)
        if idx == 0:
            return

        # STEP: accumulate r_ref to the front rnap
        self.r_ref[idx-1] += self.r_ref[idx]
        return


class DNAStrand(Entity):
    """
    This class represents the DNA strand.
    """
    def __init__(self, environment, rnap_loading_rate, include_supercoiling=True, include_busty_promoter=False,
                 rnap_loading_pattern="stochastic", promoter_shut_off_time=-1, pause_profile="flat",
                 ribo_loading_profile="stochastic", degradation_profile="exponential",
                 protein_production_off: bool = False, if_storing_supercoiling_value: bool = False,
                 implemented_t_on: float = t_on,
                 if_rnap_fall_off_from_supercoiling: bool = False, rnap_fall_off_amount: int = 5,
                 supercoiling_fall_off_upper: float = stalling_supercoiling,
                 supercoiling_fall_off_lower: float = -stalling_supercoiling):
        super().__init__(environment)
        self.length: int = length
        self.rnap_loading_rate: float = rnap_loading_rate

        # settings
        self.include_supercoiling = include_supercoiling
        self.include_busty_promoter = include_busty_promoter
        self.rnap_loading_pattern = rnap_loading_pattern
        self.pause_profile = pause_profile
        self.ribo_loading_pattern = ribo_loading_profile
        self.degradation_profile = degradation_profile
        self.protein_production_off = protein_production_off
        self.t_on = implemented_t_on
        self.if_rnap_fall_off_from_supercoiling = if_rnap_fall_off_from_supercoiling
        self.maximum_rnap_fall_off_amount = rnap_fall_off_amount
        self.rnap_fall_off_amount = 0
        self.supercoiling_fall_off_upper = supercoiling_fall_off_upper
        self.supercoiling_fall_off_lower = supercoiling_fall_off_lower

        # setup loading list.
        if_stochastic = False
        match self.rnap_loading_pattern:
            case "stochastic":
                if_stochastic = True
            case "uniform":
                if_stochastic = False
        if self.include_busty_promoter:
            self.loading_list = LoadingList(self, scaling(total_time), self.rnap_loading_rate*dt,
                                            if_stochastic=if_stochastic, if_bursty=True)
        else:
            self.loading_list = LoadingList(self, scaling(total_time), self.rnap_loading_rate*dt,
                                            if_stochastic=if_stochastic, if_bursty=False)

        # promoter_state
        self.promoter_state = False
        if promoter_shut_off_time == -1:
            self.T_stop = scaling(total_time)
        elif promoter_shut_off_time >= 0:
            self.T_stop = scaling(promoter_shut_off_time)
            self.loading_list.trim(self.T_stop)

        # Site-Specific Pausing
        self.include_site_specific_pausing = True
        match self.pause_profile:
            case "flat":
                self.include_site_specific_pausing = False
            case "OnepauseAbs":
                self.n_pause = 1
            case "TwopauseAbs":
                self.n_pause = 2

        # mRNA Degradation
        self.RNAP_LIST = RNAPList(self)

        # adaptive supercoiling
        self.r_ref = self.RNAP_LIST.r_ref
        self.flag_r_ref = self.RNAP_LIST.flag_r_ref
        self.T_open = scaling(total_time)
        self.just_loaded = False

        # storage
        self.if_storing_supercoiling_value = if_storing_supercoiling_value
        self.phi = None
        self.torq = None
        self.f_n = None
        self.velo = None
        self.stepping = None
        self.serial_number = None
        self.protein_amount = 0

    def init(self):
        self.RNAP_LIST.init()
        pass

    def step(self, time_index):

        # REASON: time to check for loading
        to_load = False

        # REASON: check if there is one loading attempt, obviously, the after T_stop, nothing should load.
        if self.loading_list.if_can_load(time_index) and time_index <= self.T_stop:
            to_load = True

        # REASON: check if there is one RNAP congesting the loading site, we just check the last rnap.
        if to_load and not self.RNAP_LIST.if_loading_site_clean():
            to_load = False

        # REASON: if it can load, then load one RNAP
        if to_load:
            if self.include_supercoiling:
                self.T_open = time_index + scaling(self.t_on)
                self.promoter_state = True
                self.just_loaded = True
                if self.RNAP_LIST.loaded > 0:
                    self.r_ref[-1] = self.RNAP_LIST.attached_rnap_list[-1].position
                    self.flag_r_ref[-1] = True

            self.RNAP_LIST.attach_rnap(initial_t=time_index, pause_profile=self.pause_profile,
                                       ribo_loading_profile=self.ribo_loading_pattern,
                                       degradation_profile=self.degradation_profile,
                                       protein_production_off=self.protein_production_off)
        # REASON: check the promoter closing resulting from the RNAP loading
        if self.just_loaded and time_index >= self.T_open and self.promoter_state:
            # REASON: switch the promoter
            self.just_loaded = False
            self.promoter_state = False
            if self.RNAP_LIST.loaded != 0:
                # REASON: we reset the reference position of the last RNAP
                self.r_ref[-1] = self.RNAP_LIST.get_rear_attached_rnap().position

        # REASON: check for permanent promoter shutoff
        if time_index >= self.T_stop and self.promoter_state:
            self.promoter_state = False
            if self.RNAP_LIST.loaded != 0:
                value = self.RNAP_LIST.attached_rnap_list[-1].position
                self.r_ref[-1] = value

        # REASON: calculate stepping
        if self.include_supercoiling:
            stepping, serial_numbers = self.supercoiling()
        else:
            stepping = []
            serial_numbers = self.RNAP_LIST.get_attached_serial_number()
            for i in range(self.RNAP_LIST.attached):
                stepping.append(v_0*dt)

        # REASON: check for site-specific pausing and set pausing.
        if self.include_site_specific_pausing and self.RNAP_LIST.attached != 0:
            for count, rnap in enumerate(self.RNAP_LIST.attached_rnap_list):
                # REASON: if passed both, then just skip
                if rnap.passed_site_1 and rnap.passed_site_2:
                    continue

                if not rnap.passed_site_1:
                    # REASON: check if the rnap is still in the pause site.
                    if rnap.passing_1 and rnap.position + 1 / scaling(pauseDuration[0]) < pauseSite[0]:
                        stepping[count] = 1 / scaling(pauseDuration[0])
                        continue

                    # REASON: passing out of the pausing site.
                    if rnap.passing_1 and rnap.position + 1 / scaling(pauseDuration[0]) >= pauseSite[0]:
                        self.RNAP_LIST.attached_rnap_list[count].passing_1 = False
                        self.RNAP_LIST.attached_rnap_list[count].passed_site_1 = True
                        continue

                    # REASON: passing into the pause_site 1
                    if rnap.position < (pauseSite[0] - 1) and ((rnap.position + stepping[count]) >= pauseSite[0] - 1):
                        self.RNAP_LIST.attached_rnap_list[count].passing_1 = True
                        self.RNAP_LIST.attached_rnap_list[count].position = pauseSite[0] - 1
                        stepping[count] = 0
                        continue

                if not rnap.passed_site_2:
                    if rnap.passing_2 and rnap.position + 1 / scaling(pauseDuration[1]) < pauseSite[1]:
                        stepping[count] = 1 / scaling(pauseDuration[1])
                        continue

                    # REASON: passing out of the pausing site.
                    if rnap.passing_2 and rnap.position + 1 / scaling(pauseDuration[1]) >= pauseSite[1]:
                        self.RNAP_LIST.attached_rnap_list[count].passing_2 = False
                        self.RNAP_LIST.attached_rnap_list[count].passed_site_2 = True
                        continue

                    # REASON: passing into the pause_site 2
                    if rnap.position < (pauseSite[1] - 1) and ((rnap.position + stepping[count]) >= pauseSite[1] - 1):
                        self.RNAP_LIST.attached_rnap_list[count].passing_2 = True
                        self.RNAP_LIST.attached_rnap_list[count].position = pauseSite[1] - 1
                        stepping[count] = 0
                        continue

        # REASON: check for hindrance and modify stepping
        for i in range(len(self.RNAP_LIST.attached_rnap_list)):
            rnap = self.RNAP_LIST.attached_rnap_list[i]
            # REASON: if this rnap is not attached, we do not need to worry about it.
            # REASON: if this rnap is the front-most attached rnap, we also do not need to care about it.
            if i == 0:
                continue
            # REASON: if this rnap is going to move past its front rnap, then a collision happens,
            #         and we need to set its position such that it sit right between its front rnap with a distance of
            #         RNAP_size.
            else:
                previous_rnap_end_position = self.RNAP_LIST.attached_rnap_list[i-1].position+stepping[i - 1]-RNAP_size
                if rnap.position + stepping[i] > previous_rnap_end_position:
                    stepping[i] = previous_rnap_end_position - rnap.position

        # REASON: now we plug in the stepping into the RNAPs and collect protein production from all RNAP.
        prot = self.RNAP_LIST.step(time_index, stepping, serial_numbers)

        # return protein production
        self.protein_amount += prot
        return prot

    def supercoiling(self):
        # STEP: setup
        positions, serial_number = self.RNAP_LIST.get_position_for_all_attached_rnap()
        r_ref, flag_r_ref = self.RNAP_LIST.get_supercoiling_ref()
        size = len(positions)
        positions = np.array(positions)

        # STEP: phi generation
        phi: ndarray = np.zeros(size+1)
        for i in range(size + 1):
            if i == 0:  # STEP: check the front-most element
                phi[i] = 0.0
            elif i == size:  # STEP: check the back-most element
                if self.promoter_state:  # STEP: check if the repressor is off. The supercoiling is diffused
                    phi[i] = 0.0
                else:  # STEP: check if the repressor is attached, full supercoiling.
                    travel_dist_last_rnap = positions[i-1] - r_ref[i - 1]
                    phi[i] = s(travel_dist_last_rnap)
            else:  # STEP: check the middle element
                travel_dist_last_rnap = positions[i - 1] - r_ref[i - 1]
                phi[i] = s(travel_dist_last_rnap - positions[i])
                # not passed the location

        # STEP: torque generation
        torq = np.zeros(size)
        n = n_dependence_cubic_3(size)

        for i in range(size):
            torq[i] = -tau_0 * n * (phi[i]-phi[i+1])

        # STEP: velocity generation
        velo = np.zeros(size)

        for i in range(size):
            # REASON: we are using this for loop, to avoid large value generated from the exponential function.
            if torq[i] > 1.5 * tau_c:
                velo[i] = 0
            elif torq[i] < -1.5 * tau_c:
                velo[i] = 2 * v_0
            else:
                velo[i] = 2 * v_0 / (1 + np.exp(2 * (torq[i] / tau_c) ** 3))
        stepping: np.ndarray = velo*dt

        # STEP: checking for RNAP fall-off due to high supercoiling
        if self.if_rnap_fall_off_from_supercoiling:
            phi, stepping, serial_number = self.rnap_fall_off(phi, stepping, serial_number)

        # STEP: relevant variables recording
        if self.if_storing_supercoiling_value:
            self.phi = phi
            self.torq = torq
            self.f_n = n
            self.velo = velo
            self.stepping = stepping
            self.serial_number = serial_number
        return stepping, serial_number

    def rnap_fall_off(self, phi, stepping, serial_number_list):
        corrected_phi = []
        corrected_stepping = []
        corrected_serial_number_list = []

        # STEP: find all the RNAPs that have high supercoiling
        high_supercoiling_rnap_serial_number_list = []
        for i in range(len(serial_number_list)):
            if if_out_of_interval(phi[i], -1, stalling_supercoiling) \
                    and self.rnap_fall_off_amount < self.maximum_rnap_fall_off_amount:
                self.rnap_fall_off_amount += 1
                high_supercoiling_rnap_serial_number_list.append(serial_number_list[i])
            else:
                corrected_phi.append(phi[i])
                corrected_stepping.append(stepping[i])
                corrected_serial_number_list.append(serial_number_list[i])

        # STEP: Remove all such rnap from the active_rnap_list
        #       Also remove the stepping of removed rnap
        for i in range(len(high_supercoiling_rnap_serial_number_list)):
            serial_n = high_supercoiling_rnap_serial_number_list[i]
            rnap = self.RNAP_LIST.get_attached_rnap(serial_n)
            self.RNAP_LIST.call_back("high_supercoiling_fall_off", rnap)

        # STEP: return the corrected stepping and serial_number_list
        return corrected_phi, corrected_stepping, corrected_serial_number_list
