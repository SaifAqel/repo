# geometry/base.py
from abc import ABC, abstractmethod
from common.units import ureg, Q_

class IGeometryStage(ABC):
    @abstractmethod
    def length(self) -> Q_:
        """Return stage length as a Quantity with length units"""
        ...

    @abstractmethod
    def hydraulic_diameter(self) -> Q_:
        """Return hydraulic diameter as a Quantity with length units"""
        ...

    @abstractmethod
    def flow_area(self) -> Q_:
        """Return flow cross-sectional area as a Quantity with area units"""
        ...

    @abstractmethod
    def heat_transfer_area(self) -> Q_:
        """Return heat transfer surface area as a Quantity with area units"""
        ...
