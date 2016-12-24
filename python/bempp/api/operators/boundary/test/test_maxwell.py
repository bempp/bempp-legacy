"""Maxwell test cases."""

from unittest import TestCase
import bempp

WAVE_NUMBER = 1


class TestMaxwell(TestCase):
    """Test cases for Maxwell operators."""

    def setUp(self):
        grid = bempp.api.shapes.regular_sphere(3)
        self._space = bempp.api.function_space(grid, "RT", 0)


if __name__ == "__main__":
    from unittest import main

    main()
