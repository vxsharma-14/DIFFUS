'''
+**************************************************************************
+**************************************************************************
+
+   FILE         readfiles.py
+
+   AUTHOR       Vishal Sharma
+
+   VERSION      1.0.0.dev1
+
+   WEBSITE      https://vxsharma-14.github.io/NAnPack/
+
+   NAnPack Learner's Edition is distributed under the MIT License.
+
+   Copyright (c) 2020 Vishal Sharma
+
+   Permission is hereby granted, free of charge, to any person
+   obtaining a copy of this software and associated documentation
+   files (the "Software"), to deal in the Software without restriction,
+   including without limitation the rights to use, copy, modify, merge,
+   publish, distribute, sublicense, and/or sell copies of the Software,
+   and to permit persons to whom the Software is furnished to do so,
+   subject to the following conditions:
+
+   The above copyright notice and this permission notice shall be
+   included in all copies or substantial portions of the Software.
+
+   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
+   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
+   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
+   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
+   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
+   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
+   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
+   SOFTWARE.
+
+   You should have received a copy of the MIT License along with
+   NAnPack Learner's Edition.
+
+**************************************************************************
+**************************************************************************
'''
#**************************************************************************
class RunConfig:
    '''This is a class to run configuration commands and read configuration
    properties from the input file with the extension ".ini".

    Attributes
    ----------

    InFileName: str, Default= './input/config.ini'.

                The string value representing the file to be read for
                simulation inputs. Use ".ini" extension files.

    '''

    def __init__(self, InFileName):
        '''The constructor for RunConfig class.

        Read input parameters for the simulation set-up from a saved
        file.

        Parameters
        ----------

        InFileName: str

                    The string value representing the file to be read for
                    simulation inputs. Use ".ini" extension files.

        '''
        import configparser
        import preprocess.checkconfig as chk

        self.File = InFileName

        print('*******************************************************')
        print('*                                                     *')
        print('*                 STARTING PRE-PROCESSING             *')
        print('*                                                     *')
        print('*******************************************************')
        print()
        print('Searching for simulation configuration file in path:')
        print(f'"{InFileName}"')
        self.config = configparser.ConfigParser()
        dataset = self.config.read(InFileName)
        if dataset:
            print('SUCCESS: Configuration file parsing.')
        else:
            raise Exception('ERROR: CONFIGURATION FILE NOT FOUND.')

        # Upon initialization -
        #       1. check - if all sections exist
        #       2. access numerical setup and and check inputs
        print('Checking whether all sections are included in config file.')
        chk.CheckSections(self.config,self.File)
        # Access all other sections and set variables.
        self.ConfigSolverSetUp()
        self.ConfigGrid()
        self.ConfigInitial()
        self.ConfigBC()
        self.ConfigConstants()
        self.ConfigSimStop()
        self.ConfigOutput()
        self.DisplayConfig()
        
#************************ FUNCTION SOLVERSETUP ****************************
    def ConfigSolverSetUp(self):
        '''Access numerical setup inputs that are specified in the
        SETUP section of configuration file.

        Call signature :

            RunConfig.SolverSetUp()
        '''
        import preprocess.checkconfig as chk
        #************ SET-UP ***************
        print('Checking numerical setup.')
        self.ExpId = self.config['SETUP']['EXPID']
        self.UnitSystem = self.config['SETUP']['UNITS_SYSTEM']
        self.Description = self.config['SETUP']['DESCRIPTION']
        self.State = self.config['SETUP']['STATE']
        self.Model = self.config['SETUP']['MODEL']
        self.Scheme = self.config['SETUP']['SCHEME']
        self.Dimension = self.config['SETUP']['DIMENSION']
        chk.CheckSetupSection(self.config,self.State,self.Model,\
                              self.Dimension,self.File)
#************************* FUNCTION GRIDGEN *******************************
    def ConfigGrid(self):
        '''Access meshing inputs that are specified in the
        DOMAIN and MESH sections of configuration file.

        Returns the grid points and grid step parameters.

        Call signature :

            RunConfig.GridGen()
        '''
        import preprocess.grid as grid
        #***************** DOMAIN SPEC *****************
        self.Length = float(self.config['DOMAIN']['LENGTH'])
        self.Height = float(self.config['DOMAIN']['HEIGHT'])
        print('Accessing domain geometry configuration: Completed')
        #**************** MESH SPEC *****************
        GridfromFile = self.config['MESH']['GRID_FROM_FILE?']
        if GridfromFile.upper() == 'YES':
            self.GridFName = self.config['MESH']['GRID_FNAME']
            if self.GridFName.lower() == 'none':
                raise Exception(f'ERROR: GRID INPUT FILE NAME.\
\nIn input file : {self.InFileName}\nIn section    : MESH\
 SPECIFICATION\nIn field      : {self.GRID_FNAME}')
            else:
                print('Functionality not available at this time')
                print('Proceeding using other inputs.')
        GridAutoCalc = self.config['MESH']['GRID_AUTO_CALC?']
        if GridAutoCalc.upper() == 'YES':
            self.dX = float(self.config['MESH']['dX'])
            self.dY = float(self.config['MESH']['dY'])
            print('Accessing meshing configuration: Completed.')
            self.iMax, self.jMax = grid.ComputeGridPoints(self.Dimension,\
                                                          self.Length,\
                                                          self.dX,\
                                                          self.Height,\
                                                          self.dY)
        elif GridAutoCalc.upper() == 'NO':
            self.iMax = int(self.config['MESH']['iMax'])
            self.jMax = int(self.config['MESH']['jMax'])
            print('Accessing meshing configuration: Completed.')
            self.dX, self.dY = grid.ComputeGridSteps(self.Dimension,\
                                                         self.Length,\
                                                         self.iMax,\
                                                         self.Height,\
                                                         self.jMax)   
#************************* FUNCTION INITIAL *******************************
    def ConfigInitial(self):
        '''Access intial condition inputs that are specified in the
        IC section of configuration file.

        Returns the dependent variables with the initial values.

        Call signature :

            RunConfig.Initial()
        '''
        import preprocess.initialize as init
        #*********** INITIAL CONDITIONS *************
        self.StartOpt = self.config['IC']['START_OPT']
        if self.StartOpt.upper() == 'RESTART':
           self.RestartFile = self.config['IC']['RESTART_FILE']
           print('Accessing initial condition settings: Completed.')
           if self.RestartFile.lower() == 'none':
               print('ERROR: SOLUTION RESTART FILE NAME.')
               print(f'In input file : {InFileName}')
               print('In section    : INITIAL CONDITIONS')
               print('In field      : RESTART_FILE')
               print('This error is generated because START = RESTART\
 option is selected without specifying the proper path to restart file.')
               ch = int(input('Do you want to proceed with:\
\n\t1. COLD START conditions, or\n\t2. Use a system default\
 RESTART filename?\nENTER 1 or 2.\n'))
               if ch == 1:
                   # initialize with zero
                   self.StartOpt = 'COLD-START'
                   self.U = init.InitialCondition(self.Dimension,\
                                                  self.iMax, self.jMax)
               elif ch == 2:
                   # use defualt file name
                   # InitFile = './output/restart.dat'
                   print('This functionality is not available at this\
 time.')
                   print('Proceeding to solve using cold-start\
 conditions.')
                   self.U  = init.InitialCondition(self.Dimension,\
                                                   self.iMax, self.jMax)

        elif self.StartOpt.upper() == 'COLD-START':
            print('Accessing initial condition settings: Completed.')
            self.U  = init.InitialCondition(self.Dimension, self.iMax,\
                                       self.jMax)

        return self.U
#**************************** FUNCTION BC *********************************
    def ConfigBC(self):
        '''Access boundary condition inputs that are specified in the
        BC section of configuration file.

        Returns the dependent variables with the boundary settings.

        Call signature :

            RunConfig.BC()
        '''
        import preprocess.boundary as bound
        #*********** BOUNDARY CONDITIONS *************
        self.BCfromFile = self.config['BC']['BC_FROM_FILE?']
        self.BCFileName = self.config['BC']['BC_FILE_NAME']
        print('Accessing boundary condition settings: Completed')
        if self.BCfromFile.upper() == 'YES':
            BC = ReadBCfromFile(self.BCFileName)
            # Call function to assign 2D BC
            self.U = bound.BC2D(self.U, BC, self.dX, self.dY)
        elif self.BCfromFile.upper() == 'NO':
            print('Follow the instructions in the documentation.')
            self.U = self.U

        return self.U
#************************* FUNCTION CONSTANTS *****************************
    def ConfigConstants(self):
        '''Access constant inputs that are specified in the
        CONST section of configuration file.

        Call signature :

            RunConfig.Constants()
        '''
        #*********** CONSTANT COEFFICIENTS ***********
        self.CFL = float(self.config['CONST']['CFL'])
        self.conv = float(self.config['CONST']['CONV'])
        self.diff = float(self.config['CONST']['DIFF'])
        print('Accessing constant data: Completed.')
#*********************** FUNCTION CONFIGSIMSTOP ***************************
    def ConfigSimStop(self):
        '''Access simulation stop setting inputs that are specified in the
        STOP section of configuration file.

        Call signature :

            RunConfig.ConfigSimStop()
        '''
        import preprocess.time as time
        #*********** SIM STOP SETTINGS ***********
        self.totTime = float(self.config['STOP']['SIM_TIME'])
        if self.State.upper() == 'STEADY':
            self.ConvCrit = float(self.config['STOP']['CONV_CRIT'])
        elif self.State.upper() == 'TRANSIENT':
            self.ConvCrit = -0.01
        nMax = int(self.config['STOP']['nMAX'])
        # Execute this block for:
        # diffusion eq., first-order wave eq. and Burgers eq.
        if not self.Model.upper() == 'POISSONS':
            self.dT = time.CalcTimeStep(self.CFL, self.diff, self.conv,\
                                        self.dX, self.dY,\
                                        self.Dimension, self.Model)
            self.nMax = time.CalcMaxSteps(self.State, nMax, self.dT,\
                                          self.totTime)
        # Execute this block for Poissons eq.
        elif self.Model.upper() == 'POISSONS':
            self.nMax = nMax
            
        print('Accessing simulation stop settings: Completed.')
#*********************** FUNCTION CONFIGOUTPUT ****************************
    def ConfigOutput(self):
        '''Access output configurations that are specified in the
        OUTPUT section of configuration file.

        Call signature :

            RunConfig.ConfigOutput()
        '''
        #************* OUTPUT INFORMATION *************        
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
                nodeX =  self.config['OUTPUT']['X']
                self.nodes = [float(node) for node in nodeX.split(',')]
                self.PrintNodesDir = 'X'
            except:
                nodeY =  self.config['OUTPUT']['Y']
                self.nodes = [float(node) for node in nodeY.split(',')]
                self.PrintNodesDir = 'Y'
            self.Out1DFName = self.config['OUTPUT']['SAVE1D_FILENAME']
        print('Accessing settings for storing outputs: Completed.')
#*********************** FUNCTION DISPLAYCONFIG ***************************
    def DisplayConfig(self):
        '''Display saved configuration to the user for verification.

        Call signature :

            RunConfig.DisplayConfig()
        '''
        #************* PRINT CONFIGURATIONS *************
        print()
        print('***********************************************************')
        print(f'CASE DESCRIPTION                {self.Description}')
        print(f'SOLVER STATE                    {self.State}')
        print(f'MODEL EQUATION                  {self.Model}')
        print(f'DOMAIN DIMENSION                {self.Dimension}')
        print(f'    LENGTH                      {self.Length}')
        if self.Dimension.upper() == '2D':
            print(f'    HEIGHT                      {self.Height}')
        print('GRID STEP SIZE')
        print(f'    dX                          {self.dX:5.3f}')
        if self.Dimension.upper() == '2D':
            print(f'    dY                          {self.dY:5.3f}')
        if not self.Model == 'POISSONS':
            print(f'TIME STEP                   {self.dT:5.3f}')
        print('GRID POINTS')
        print(f'    along X                     {self.iMax}')
        if self.Dimension.upper() == '2D':
            print(f'    along Y                     {self.jMax}')
        if self.Model.upper() == 'DIFFUSION':
            print(f'DIFFUSION CONST.                {self.diff:6.4e}')
            print(f'DIFFUSION NUMBER                {self.CFL}')
        elif self.Model.upper() == 'FO_WAVE':
            print(f'CONVECTION CONST.               {self.conv}')
            print(f'COURANT NUMBER                  {self.CFL}')
        elif self.Model.upper() == 'BURGERS':
            print(f'COURANT NUMBER                  {self.CFL}')
        if self.State.upper() == 'STEADY':
            print(f'CONVERGENCE CRIT                {self.ConvCrit}')
            print(f'MAXIMUM ITERATIONS')
            print(f'IF CONVERGENCE NOT OBSERVED     {self.nMax}')
        elif self.State.upper() == 'TRANSIENT':
            print(f'TOTAL SIMULATION TIME           {self.totTime}')
            print(f'NUMBER OF TIME STEPS            {self.nMax}')
        if self.BCfromFile.upper() == 'YES':
            print(f'BC COFIGURATION FILE            "{self.BCFileName}"')
        print(f'START CONDITION                 {self.StartOpt}')
        print('***********************************************************')
        print('SUCEESS: Configuration completed.')
        print()
#**************************************************************************
def ReadBCfromFile(BCFileName):
    '''Read user specified boundary conditions from the file

    Call signature;

        ReadBCfromFile(BCFileName)

    Parameters
    __________

    InFileName: str

                The string value representing the file to be read for
                simulation inputs.
    '''
    import configparser
    import preprocess.checkconfig as chk

    print('***********************************************************')
    print(f'READING BOUNDARY CONDITIONS CONFIGURATION FROM FILE:')
    print(f'"{BCFileName}"')
    print('***********************************************************')

    config = configparser.ConfigParser()
    dataset = config.read(BCFileName)

    if dataset:
        print('SUCCESS: Boundary conditions configuration file parsing.')
    else:
        raise Exception('BC config file not found :\n{BCFileName}.')

    sections = ['INLET', 'WALL', 'FAR-FIELD', 'OUTLET']
    #*********** CHECK IF ALL SECTIONS EXIST *********
    print('Checking whether all sections are included in config file.')
    chk.CheckBCSections(config,BCFileName)

    #************ INLET BC ***************
    InletAxis1 = config['INLET']['AXIS_1']
    Ain1 = float(config['INLET']['A1'])
    Bin1 = float(config['INLET']['B1'])
    InletAxis2 = config['INLET']['AXIS_2']
    if InletAxis2.lower() != 'none':
        Ain2 = float(config['INLET']['A2'])
        Bin2 = float(config['INLET']['B2'])
    else:
        Ain2 = 'none'
        Bin2 = 'none'
    InBCType = config['INLET']['BC_TYPE']
    Uin = float(config['INLET']['U'])
    print('Accessing boundary conditions at INLET: Completed.')

    #************ WALL BC **************
    WallAxis1 = config['WALL']['AXIS_1']
    Aw1 = float(config['WALL']['A1'])
    Bw1 = float(config['WALL']['B1'])
    WallAxis2 = config['WALL']['AXIS_2']
    if WallAxis2.lower() != 'none':
        Aw2 = float(config['WALL']['A2'])
        Bw2 = float(config['WALL']['B2'])
    else:
        Aw2 = 'none'
        Bw2 = 'none'
    WallBCType = config['WALL']['BC_TYPE']
    Uw = float(config['WALL']['U'])
    print('Accessing boundary conditions at WALL: Completed.')

    #************* FAR-FIELD BC ***************
    FarfldAxis1 = config['FAR-FIELD']['AXIS_1']
    Aff1 = float(config['FAR-FIELD']['A1'])
    Bff1 = float(config['FAR-FIELD']['B1'])
    FarfldAxis2 = config['FAR-FIELD']['AXIS_2']
    if FarfldAxis2.lower() != 'none':
        Aff2 = float(config['FAR-FIELD']['A2'])
        Bff2 = float(config['FAR-FIELD']['B2'])
    else:
        Aff2 = 'none'
        Bff2 = 'none'
    FarfldBCType = config['FAR-FIELD']['BC_TYPE']
    Uff = float(config['FAR-FIELD']['U'])
    print('Accessing boundary conditions at FAR-FIELD: Completed.')

    #*********** OUTLET BC **************
    OutletAxis1 = config['OUTLET']['AXIS_1']
    Ao1 = float(config['OUTLET']['A1'])
    Bo1 = float(config['OUTLET']['B1'])
    OutletAxis2 = config['OUTLET']['AXIS_2']
    if OutletAxis2.lower() != 'none':
        Ao2 = float(config['OUTLET']['A2'])
        Bo2 = float(config['OUTLET']['B2'])
    else:
        Ao2 = 'none'
        Bo2 = 'none'
    OutBCType = config['OUTLET']['BC_TYPE']
    Uo = float(config['OUTLET']['U'])
    print('Accessing boundary conditions at OUTLET: Completed.')

    print('Reading and acessing BC from config file: Completed.')
    Inlet = [InletAxis1, Ain1, Bin1, InletAxis2, Ain2, Bin2, InBCType, Uin]
    Wall = [WallAxis1, Aw1, Bw1, WallAxis2, Aw2, Bw2, WallBCType, Uw]
    Farfld = [FarfldAxis1, Aff1, Bff1, FarfldAxis2, Aff2, Bff2,\
              FarfldBCType, Uff]
    Outlet = [OutletAxis1, Ao1, Bo1, OutletAxis2, Ao2, Bo2, OutBCType, Uo]

    print('***********************************************************')

    return Inlet, Wall, Farfld, Outlet 

#**************************************************************************
'''if __name__ == '__main__':
    import sys
    ReadConfig(sys.argv[1])
'''