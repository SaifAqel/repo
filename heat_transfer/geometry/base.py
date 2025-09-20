# geometry/base.py
from abc import ABC, abstractmethod
class IGeometryStage(ABC):
    @abstractmethod
    def length(self) -> float:
        ...
    @abstractmethod
    def hydraulic_diameter(self) -> float:
        ...
    @abstractmethod
    def flow_area(self) -> float:
        ...
    @abstractmethod
    def heat_transfer_area(self) -> float:
        ...
