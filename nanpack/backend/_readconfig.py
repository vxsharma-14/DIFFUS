"""."""
import configparser
from .exceptions import InputFileError


class _ReadConfig:
    """Class to read inputs from configuration file.

    This class is depenedent on the configparser package to read the
    inputs from the file and therefore, for compability,
    always use the most recent version of the config file provided
    with the package.

    Attributes
    ----------
    InFileName: str, Default= "./input/config.ini".

        The string value representing the file to be read for simulation
        inputs. Allowed file extensions: .ini
    """

    def __init__(self, InFileName):
        """Class constructor for the _ReadConfig class.

        Parameters
        ----------
        InFileName: str

            Path to configuration file.
        """
        self.File = InFileName
        self.config = configparser.ConfigParser()
        dataset = self.config.read(InFileName)
        if not dataset:
            raise InputFileError("FileNotFound", InFileName)

        self.ExpId = self.config['SETUP']['EXPID']
        self.UnitSystem = self.config['SETUP']['UNITS_SYSTEM']
        self.Description = self.config['SETUP']['DESCRIPTION']
        self.State = self.config['SETUP']['STATE']
        self.Model = self.config['SETUP']['MODEL']
        self.Scheme = self.config['SETUP']['SCHEME']
        self.Dimension = self.config['SETUP']['DIMENSION']
        # ***************** DOMAIN SPEC *****************
        self.Length = float(self.config['DOMAIN']['LENGTH'])
        self.Height = float(self.config['DOMAIN']['HEIGHT'])
        # **************** MESH SPEC *****************
        self.GridfromFile = self.config['MESH']['GRID_FROM_FILE?']
        if self.GridfromFile.upper() == 'YES':
            self.GridFName = self.config['MESH']['GRID_FNAME']
            if self.GridFName.lower() == 'none':
                raise InputFileError("FileNotFound", self.GridFName)
            else:
                print("Functionality not available at this time")
                print("Proceeding using other inputs.")
        self.GridAutoCalc = self.config['MESH']['GRID_AUTO_CALC?']
        if self.GridAutoCalc.upper() == 'YES':
            self.dX = float(self.config['MESH']['dX'])
            self.dY = float(self.config['MESH']['dY'])
        elif self.GridAutoCalc.upper() == 'NO':
            self.iMax = int(self.config['MESH']['iMax'])
            self.jMax = int(self.config['MESH']['jMax'])

        self.StartOpt = self.config['IC']['START_OPT']
        if self.StartOpt.upper() == 'RESTART':
            print("Functionality does not exists in this version.")
            pass
        if 1 == 0:
            self.RestartFile = self.config['IC']['RESTART_FILE']
            print("Accessing initial condition settings: Completed.")
            if self.RestartFile.lower() == "none":
                raise InputFileError("FileNotFound", self.RestartFile)
                ch = int(input("Do you want to proceed with:\
\n\t1. COLD START conditions, or\n\t2. Use a system default\
 RESTART filename?\nENTER 1 or 2.\n"))
                if ch == 1:
                    # Initialize with zero
                    self.StartOpt = "COLD-START"
                elif ch == 2:
                    # Use defualt file name
                    # InitFile = './output/restart.dat'
                    print("This functionality is not available at this\
 time.")
                    print("Proceeding to solve using cold-start\
 conditions.")
        # *********** BOUNDARY CONDITIONS *************
        self.BCfromFile = self.config['BC']['BC_FROM_FILE?']
        self.BCFileName = self.config['BC']['BC_FILE_NAME']
        # *********** CONSTANT COEFFICIENTS ***********
        self.CFL = float(self.config['CONST']['CFL'])
        self.conv = float(self.config['CONST']['CONV'])
        self.diff = float(self.config['CONST']['DIFF'])
        # *********** SIM STOP SETTINGS ***********
        self.totTime = float(self.config['STOP']['SIM_TIME'])
        if self.State.upper() == 'STEADY':
            self.ConvCrit = float(self.config['STOP']['CONV_CRIT'])
        elif self.State.upper() == 'TRANSIENT':
            self.ConvCrit = -0.01
        self.nMax = int(self.config['STOP']['nMAX'])
        # ************* OUTPUT INFORMATION *************
        self.HistFileName = self.config['OUTPUT']['HIST_FILE_NAME']
        self.RestartFile = self.config['OUTPUT']['RESTART_FNAME']
        self.OutFileName = self.config['OUTPUT']['RESULT_FNAME']
        self.nWrite = int(self.config['OUTPUT']['WRITE_EVERY'])
        self.nDisplay = int(self.config['OUTPUT']['DISPLAY_EVERY'])
        self.SaveforAnim = self.config['OUTPUT']['SAVE_FOR_ANIM?']
        if self.SaveforAnim.upper() == 'YES':
            self.nAnime = int(self.config['OUTPUT']['SAVE_EVERY'])
        self.Save1DOut = self.config['OUTPUT']['SAVE_1D_OUTPUT?']
        if self.Save1DOut.upper() == 'YES':
            nodeX = None
            nodeY = None
            try:
                nodeX = self.config['OUTPUT']['X']
                self.nodes = [float(node) for node in nodeX.split(',')]
                self.PrintNodesDir = 'X'
            except:
                nodeY = self.config['OUTPUT']['Y']
                self.nodes = [float(node) for node in nodeY.split(',')]
                self.PrintNodesDir = 'Y'
            self.Out1DFName = self.config['OUTPUT']['SAVE1D_FILENAME']
