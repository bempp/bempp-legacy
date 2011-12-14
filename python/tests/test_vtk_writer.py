#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import the bempp module
import sys
sys.path.append("..")
import bempp
import numpy as np
import os
import shutil

import pytest

class TestVtkWriter:
    def setup_method(self, method):
        self.grid = bempp.GridFactory.createStructuredGrid("triangular", (0., 0.), (1., 2.), (4, 5))
        self.view = self.grid.leafView()
        self.vtkWriter = self.view.vtkWriter()


    def test_addCellData_throws_for_input_of_incorrect_dimensions(self):
        data = np.zeros((3, 4))
        pytest.raises((RuntimeError, ValueError, TypeError),
            "self.vtkWriter.addCellData(data, 'title')")

    def test_addCellData_does_not_throw_for_input_of_correct_dimensions(self):
        data = np.zeros((self.view.entityCount(0), 4))
        pytest.raises((RuntimeError, ValueError, TypeError),
            "self.vtkWriter.addCellData(data, 'title')")

    def test_addCellData_saves_a_ref_to_python_data(self):
        data = np.zeros((4, self.view.entityCount(0)))
        self.vtkWriter.addCellData(data, 'title')
        assert len(self.vtkWriter._python_data) == 1

    def test_addCellData_saves_a_ref_to_correct_python_data(self):
        data = np.zeros((4, self.view.entityCount(0)))
        self.vtkWriter.addCellData(data, 'title')
        assert np.all(self.vtkWriter._python_data[0] == data)

    def test_addCellData_saves_a_ref_to_c_data(self):
        data = np.zeros((4, self.view.entityCount(0)))
        self.vtkWriter.addCellData(data, 'title')
        assert len(self.vtkWriter._c_data) == 1


    def test_addVertexData_throws_for_input_of_incorrect_dimensions(self):
        data = np.zeros((3, 4))
        pytest.raises((RuntimeError, ValueError, TypeError),
            "self.vtkWriter.addVertexData(data, 'title')")

    def test_addVertexData_does_not_throw_for_input_of_correct_dimensions(self):
        data = np.zeros((4, self.view.entityCount(2)))
        vtkWriter = self.vtkWriter.addVertexData(data, 'title') # should not throw

    def test_addVertexData_saves_a_ref_to_python_data(self):
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        assert len(self.vtkWriter._python_data) == 1

    def test_addVertexData_saves_a_ref_to_correct_python_data(self):
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        assert np.all(self.vtkWriter._python_data[0] == data)

    def test_addVertexData_saves_a_ref_to_c_data(self):
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        assert len(self.vtkWriter._c_data) == 1


    def test_clear_clears_refs_to_python_data(self):
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        self.vtkWriter.clear()     
        assert ((not hasattr(self.vtkWriter, "_python_data")) 
            or self.vtkWriter._python_data is None)

    def test_clear_clears_refs_to_c_data(self):
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        self.vtkWriter.clear()     
        assert ((not hasattr(self.vtkWriter, "_c_data")) 
            or self.vtkWriter._c_data is None)

    @pytest.mark.parametrize("type", ("ascii", "base64", "appendedraw", "appendedbase64"))
    def test_write_accepts_output_type(self, type):
        if os.path.exists("output.vtu"):
            os.remove("output.vtu")
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        self.vtkWriter.write("output", type)
        assert os.path.exists("output.vtu")
        os.remove("output.vtu")


    @pytest.mark.parametrize("type", ("ascii", "base64", "appendedraw", "appendedbase64"))
    def test_pwrite_accepts_output_type(self, type):
        if os.path.exists("subdir"):
            shutil.rmtree("subdir")
        os.mkdir("subdir")
        os.mkdir("subdir/nested_subdir")
        data = np.zeros((4, self.view.entityCount(2)))
        self.vtkWriter.addVertexData(data, 'title')
        self.vtkWriter.pwrite("output", "subdir", "nested_subdir", type)
        assert os.path.exists("subdir/s0001:output.pvtu")
        assert os.path.exists("subdir/nested_subdir/s0001:p0000:output.vtu")
        shutil.rmtree("subdir")