class Measurement:
    """Class modeling a measurment on the power grid
    categories can be mapped from 0 to 5 as:
        0: Active power Injection at bar(len(buses)=1)
        1: Reactive power Injection at bar(len(buses)=1)
        2: Active power Flow between bars(len(buses)=2)
        3: Reactive power Flow between bars(len(buses)=2)
        4: Voltage Level at bar(len(buses)=1)
        5: Current Flow between bars(len(buses)=2)"""
    
    def __init__(self, category: int, bus: tuple, value: float):
        
        self.category = category
        self.bus = bus
        self.value = value
    
class StateData:
    """"Class for modelling and storing data from each state of the system state vector
    categories are defined as such:
        0: angle variable
        1: voltage variable"""
        
    def __init__(self, category: int, index: int, bus):
        self.category = category
        self.index = index
        self.bus = bus
        self.value = []
        
    def addValue(self, value):
        self.value.append(value)