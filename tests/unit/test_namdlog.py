import unittest
import pytest
from pestifer.namdlog import NAMDLog,getinfo

class TestNAMDLog(unittest.TestCase):
    def test_namdlog_init(self):
        L=NAMDLog('test_namdlog.testlog')
        self.assertEqual(L.numlines,291327)
        self.assertEqual(len(L.groups['Info:']),10227)

    def test_namdlog_edata(self):
        L=NAMDLog('test_namdlog.testlog').energy()
        self.assertEqual(L.edata.shape[0],100001)
        etitles='TS           BOND          ANGLE          DIHED          IMPRP               ELECT            VDW       BOUNDARY           MISC        KINETIC               TOTAL           TEMP      POTENTIAL         TOTAL3        TEMPAVG            PRESSURE      GPRESSURE         VOLUME       PRESSAVG      GPRESSAVG'.split()
        self.assertTrue(all([L.edata.columns[i]==etitles[i] for i in range(len(etitles))]))

    def test_namdlog_success(self):
        L=NAMDLog('test_namdlog.testlog').energy()
        self.assertTrue(L.success())

    def test_namdlog_getinfo_actual(self):
        L=NAMDLog('test_namdlog.testlog')
        self.assertTrue(L.info['TOTAL MASS']=='1.51197e+06')

    def test_namdlog_getinfo_mock(self):
        line='Info: MARGIN                 0'
        val=getinfo('MARGIN',line)
        self.assertEqual(val,'0')
        line='Info: 245567 ATOMS'
        val=getinfo('ATOMS',line)
        self.assertEqual(val,'245567')
        line='Info: 4 ATOMS IN LARGEST HYDROGEN GROUP'
        val=getinfo('ATOMS IN LARGEST HYDROGEN GROUP',line)
        self.assertEqual(val,'4')
        line='Info: ATOM DENSITY = 0.100831 atoms/A^3'
        val=getinfo('ATOM DENSITY',line)
        self.assertEqual(val,'0.100831')
