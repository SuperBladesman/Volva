"""
Standard Reference Materials (SRM) database for XRD analysis.

Contains crystallographic data and expected peak positions for common
XRD standard reference materials used for instrumental broadening correction.
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from ..core.constants import WAVELENGTHS


@dataclass
class CrystalStructure:
    """Crystallographic structure information for SRM materials."""
    name: str
    formula: str
    crystal_system: str
    space_group: str
    lattice_parameters: Dict[str, float]  # a, b, c, alpha, beta, gamma
    density: float  # g/cm³
    reference: str


@dataclass
class ReflectionData:
    """Individual reflection data for SRM peaks."""
    hkl: Tuple[int, int, int]
    d_spacing: float  # Å
    relative_intensity: float  # 0-100
    multiplicity: int
    structure_factor_squared: float


class StandardReferenceMaterial:
    """Complete SRM data including structure and expected reflections."""

    def __init__(self, structure: CrystalStructure, reflections: List[ReflectionData]):
        self.structure = structure
        self.reflections = reflections
        self._peak_positions = {}  # Cache for different wavelengths

    def get_peak_positions(self, wavelength: float, min_intensity: float = 1.0,
                           max_two_theta: float = 150.0) -> np.ndarray:
        """
        Calculate 2θ positions for reflections above minimum intensity.

        Args:
            wavelength: X-ray wavelength in Å
            min_intensity: Minimum relative intensity to include
            max_two_theta: Maximum 2θ angle to consider

        Returns:
            Array of 2θ positions in degrees
        """
        cache_key = (wavelength, min_intensity, max_two_theta)
        if cache_key in self._peak_positions:
            return self._peak_positions[cache_key]

        positions = []
        for refl in self.reflections:
            if refl.relative_intensity >= min_intensity:
                # Bragg's law: nλ = 2d sinθ
                sin_theta = wavelength / (2 * refl.d_spacing)
                if sin_theta <= 1.0:  # Valid reflection
                    theta_rad = np.arcsin(sin_theta)
                    two_theta = 2 * np.degrees(theta_rad)
                    if two_theta <= max_two_theta:
                        positions.append(two_theta)

        result = np.array(sorted(positions))
        self._peak_positions[cache_key] = result
        return result

    def get_reflection_info(self, two_theta: float, wavelength: float,
                            tolerance: float = 0.1) -> Optional[ReflectionData]:
        """
        Get reflection information for a given 2θ position.

        Args:
            two_theta: 2θ position in degrees
            wavelength: X-ray wavelength in Å
            tolerance: Angular tolerance in degrees

        Returns:
            ReflectionData if found within tolerance, None otherwise
        """
        for refl in self.reflections:
            sin_theta = wavelength / (2 * refl.d_spacing)
            if sin_theta <= 1.0:
                expected_2theta = 2 * np.degrees(np.arcsin(sin_theta))
                if abs(two_theta - expected_2theta) <= tolerance:
                    return refl
        return None


# Define standard reference materials
def create_lab6_srm() -> StandardReferenceMaterial:
    """Create LaB₆ (NIST SRM 660a) standard reference material."""
    structure = CrystalStructure(
        name="Lanthanum Hexaboride",
        formula="LaB₆",
        crystal_system="cubic",
        space_group="Pm-3m",
        lattice_parameters={"a": 4.15695, "b": 4.15695, "c": 4.15695,
                            "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        density=4.72,
        reference="NIST SRM 660a"
    )

    # LaB₆ reflection data (major peaks)
    reflections = [
        ReflectionData((1, 0, 0), 4.157, 100.0, 6, 1000.0),
        ReflectionData((1, 1, 0), 2.939, 52.3, 12, 520.0),
        ReflectionData((1, 1, 1), 2.400, 26.8, 8, 268.0),
        ReflectionData((2, 0, 0), 2.078, 20.1, 6, 201.0),
        ReflectionData((2, 1, 0), 1.859, 22.5, 24, 225.0),
        ReflectionData((2, 1, 1), 1.699, 15.8, 24, 158.0),
        ReflectionData((2, 2, 0), 1.470, 14.2, 12, 142.0),
        ReflectionData((3, 0, 0), 1.386, 8.9, 6, 89.0),
        ReflectionData((3, 1, 0), 1.316, 12.7, 24, 127.0),
        ReflectionData((3, 1, 1), 1.254, 11.3, 24, 113.0),
        ReflectionData((2, 2, 2), 1.200, 5.2, 8, 52.0),
    ]

    return StandardReferenceMaterial(structure, reflections)


def create_silicon_srm() -> StandardReferenceMaterial:
    """Create Silicon (NIST SRM 640c) standard reference material."""
    structure = CrystalStructure(
        name="Silicon",
        formula="Si",
        crystal_system="cubic",
        space_group="Fd-3m",
        lattice_parameters={"a": 5.431123, "b": 5.431123, "c": 5.431123,
                            "alpha": 90.0, "beta": 90.0, "gamma": 90.0},
        density=2.329,
        reference="NIST SRM 640c"
    )

    # Silicon reflection data (diamond structure, only even h+k+l)
    reflections = [
        ReflectionData((1, 1, 1), 3.135, 100.0, 8, 1000.0),
        ReflectionData((2, 2, 0), 1.920, 55.0, 12, 550.0),
        ReflectionData((3, 1, 1), 1.637, 30.0, 24, 300.0),
        ReflectionData((4, 0, 0), 1.358, 6.0, 6, 60.0),
        ReflectionData((3, 3, 1), 1.246, 11.0, 24, 110.0),
        ReflectionData((4, 2, 2), 1.109, 12.0, 24, 120.0),
        ReflectionData((3, 3, 3), 1.045, 3.0, 8, 30.0),
        ReflectionData((5, 1, 1), 1.045, 3.0, 24, 30.0),
        ReflectionData((4, 4, 0), 0.960, 9.0, 12, 90.0),
        ReflectionData((5, 3, 1), 0.918, 8.0, 48, 80.0),
    ]

    return StandardReferenceMaterial(structure, reflections)


# SRM Registry
STANDARD_MATERIALS = {
    "LaB6": create_lab6_srm(),
    "Lab6": create_lab6_srm(),  # Alternative name
    "Silicon": create_silicon_srm(),
    "Si": create_silicon_srm(),  # Alternative name
}


def get_srm(name: str) -> Optional[StandardReferenceMaterial]:
    """
    Retrieve a standard reference material by name.

    Args:
        name: Name of the SRM (case-insensitive)

    Returns:
        StandardReferenceMaterial instance or None if not found
    """
    return STANDARD_MATERIALS.get(name.replace(" ", ""))


def list_available_srms() -> List[str]:
    """Get list of available SRM names."""
    return ["LaB6", "Silicon"]


def recommend_srm_for_range(two_theta_range: Tuple[float, float],
                            wavelength: float) -> List[str]:
    """
    Recommend suitable SRMs based on the 2θ measurement range.

    Args:
        two_theta_range: (min_2theta, max_2theta) in degrees
        wavelength: X-ray wavelength in Å

    Returns:
        List of recommended SRM names
    """
    min_2theta, max_2theta = two_theta_range
    recommendations = []

    for name in list_available_srms():
        srm = get_srm(name)
        if srm:
            peaks = srm.get_peak_positions(wavelength, min_intensity=5.0,
                                           max_two_theta=max_2theta)
            # Check if SRM has good coverage in the range
            peaks_in_range = peaks[(peaks >= min_2theta) & (peaks <= max_2theta)]
            if len(peaks_in_range) >= 3:  # At least 3 peaks for good fitting
                recommendations.append(name)

    return recommendations