#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2013.04.16.
Under GPL licence.

Purpose:
========
Generate nmrpipe processing scripts for 1D, 2D or 3 dimensional NMR data

Requires:
=========
* NMR data in varian format
* nmrpipe processing program

"""

import os
import re
import sys

class ProcparData():
    """
    Read a varian procpar file and load the data into a variables

    Example:

        MyProcpar = ProcparData('data/','procpar')
        print MyProcpar.parameter('ni')
    """
    #####################################################
    def __init__(self, Path='', FileName='procpar'):
        """
        Parameters:

        * Path      Path of the procpar file
        * FileName  Filename of the procpar file, default = 'procpar'
        """
        # Open and read the procpar file
        try:
            File_handler = open(Path+FileName,'r')
        except IOError:
            print ''.join(('\n-----------\nError opening ', Path, FileName,
                           '! Please check it!\n-----------\n'))
            exit()
        Lines = File_handler.readlines()
        File_handler.close()
        # Data stored in __parameter__
        self.__parameter__ = {}
        # Go through the lines
        i = 0
        while (i < len(Lines)):
            # The 'key' is the first line first value
            key = Lines[i].split()[0]
            while key == '0':
                key = Lines[i].split()[0]
                i += 1
#            print '\n',key
            if (key == 'saveglobal_') or (key == 'dgs') or (key == 'dg2'):
                # The 'values' are the next lines till a "0" comes
                i += 1
                first_row = 1
                while (i < len(Lines)) and (Lines[i].strip() != '0'):
                    if first_row == 1:
                        # The first value in the first row is the numer of
                        # values, do not need to store
                        value = Lines[i].strip().split()[1:]
                        first_row = 0
                    else:
                        value.append(Lines[i].strip())
                    i += 1
                # The final line for each parameter is the "0"
                if (i < len(Lines)):
                    i += 1
            else:
                value = []
                value.append(Lines[i+1].strip()[Lines[i+1].index(' ')+1:])

                i += 3
            # Store the data
#            print key, value
            self.__parameter__[key] = value
            key = ''
        return None
    #####################################################
    def parameter(self, ParameterName):
        """
        Returns the value of a parameter from the procpar file

        Parameters:

        * ParameterName = The name of a parameter, like: ni, nt, sw2,...
        """
        if ParameterName in self.__parameter__:
            returnvalue = self.__parameter__[ParameterName]
        else:
            print '\n------------\nNo parameter like "'+ParameterName+'" in the procpar file\n------------\n'
            #exit()
            returnvalue = None
        return returnvalue
################################################################################


class Convert_HSQC():
    """
    Note:
        * Last value must be the proton phase correction
    """
    def __init__(self, argumentlist, FidFileName = 'fid', Path = './'):
        #
        if 'help' in argumentlist:
            print ('\n-------------------------------------------------\n'
                   'Usage: hsqc.com [temp] [noplot] <proton p0 phase>\n'
                   '-------------------------------------------------\n'
                   'temp      = use "temp" flag to set a different measurement temperature\n'
                   'noplot    = use "noplot" flag not to see nmrDraw and Sparky plots\n'
                   'nocleanup = keeping all files\n'
                   'extract   = proton dimension extraction, next parameter must be [6.8,10.0]\n'
                   'p0 phase  = proton phase correction value must be the last parameter\n')
            exit()
        #
        self.Path            = Path
        self.FidFileName     = FidFileName
        self.ConvertFileName = 'convert_nmr.com'
        #
        self.__temporary_folder = 'data'
        #
        self.PropcarInformation = ProcparData(Path = Path)
        #
        #
        if len(argumentlist) > 1:
            protonphase = argumentlist[-1]
        else:
            protonphase = '0.0'
        self.temp = self.Info('temp')[0]
        if 'temp' in argumentlist:
            self.SetTemperature()
        ########
        if self.Info('ni') and self.Info('ni2'):
            self._2D = self.Info('ni')[0] != '1' or self.Info('ni2')[0] != '1'
        else:
            self._2D = False
        ######
        if 'extract' in argumentlist:
            ex = eval(argumentlist[argumentlist.index('extract') + 1])
        else:
            if self._2D:
                #TODO What if it is a caron hsqc?
                #TODO The second dimension question is not solved yet!!!
                #ex = [-1.0, 4.0]
                ex = [15.0, 5.0]
                #ex = None
            else:
                # If it is a 1D then only extract if needed!
                ex = None
        #
        #
        self.CreateConverFile(userphase = protonphase, SecondDimension='N', Fastprocess='fast' in argumentlist, Extract=ex,
                              Trosy_experiment = 'trosy' in self.Info('seqfil')[0])
        #
        self.RunConvertFile(not 'noplot' in argumentlist, not 'noplot' in argumentlist, 'nocleanup' in argumentlist)
        #
        self.ByeBye()
        return None
    ###################
    def Info(self, paramatername):
        return self.PropcarInformation.parameter(paramatername)
    ###################
    def Get_H2O_chemical_shift(self, temperature):
        """
        Calculate the carrier frequency based on the temperature
        """
        return 5.01165675545 + temperature * -0.00955018024414
    ###################
    def Get_Current_Dir(self):
        """
        """
        return os.getcwd().split('/')[-1]
    ###################
    def SetTemperature(self):
        self.temp = ''
        print '\n_______T E M P E R A T U R E :_______________'
        self.temp = raw_input('Please enter the correct temperature: ')
        print 'The measurement temperature is set to ' + self.temp + ' C => XCAR: {0:8.6f}'.format(self.Get_H2O_chemical_shift(float(self.temp)))
        self.temp = float(self.temp)
        return None
    ###################
    def Get_Previous_Phase_values(self):
        Xp0 = 150.0
        Xp1 = 0.0
        Yp0 = 0.0
        Yp1 = 0.0
        if os.path.exists(self.Path + self.ConvertFileName):
            fil = open(self.Path + self.ConvertFileName,'r')
            lines = fil.readlines()
            fil.close()
            first = True
            for line in lines:
                if re.search('^\| nmrPipe -fn PS',line):
                    if first:
                        try:
                            Xp0 = float(line.split()[5])
                        except IndexError:
                            pass
                        except ValueError:
                            pass
                        try:
                            Xp1 = float(line.split()[7])
                        except IndexError:
                            pass
                        except ValueError:
                            pass
                        first = False
                    else:
                        try:
                            Yp0 = float(line.split()[5])
                        except IndexError:
                            pass
                        except ValueError:
                            pass
                        try:
                            Yp1 = float(line.split()[7])
                        except IndexError:
                            pass
                        except ValueError:
                            pass
        return [Xp0,Xp1,Yp0,Yp1]
    ###################
    def Get_carrier_in_PPM(self, Proton_carrier_PPM, Frequency_H, Frequency_X, Nucleus_type_X):
        """
        Returns the carrier in ppm

        1H, 13C and 15N chemical shift referencing in biomolecular NMR.
        Wishart DS, Bigam CG, Yao J, Abildgaard F, Dyson HJ, Oldfield E, Markley JL, Sykes BD.
        J Biomol NMR. 1995 Sep;6(2):135-40. PMID: 8589602

        http://www.bmrb.wisc.edu/ref_info/cshift.html
        """
        RATIOS = {'H':1.0, 'N':0.101329118, 'C':0.251449530}
        #
        return Frequency_X * (1E+6 + Proton_carrier_PPM) / (RATIOS[Nucleus_type_X] * Frequency_H) - 1E+6
    ###################
    def script_for_regular_1D(self, parameters):
        """
        """
        script = ('#!/bin/csh\n'
                  '\n'
                  'var2pipe -in {_fidfile_} -noaswap \\\n'
                  '     -xN    {_xN____:>12s}        \\\n'
                  '     -xT    {_xT____:>12s}        \\\n'
                  '     -xMODE {_xMODE_:>12s}        \\\n'
                  '     -xSW   {_sw____:>12s}        \\\n'
                  '     -xOBS  {_xobs__:>12s}        \\\n'
                  '     -xCAR  {_XCAR__:>12s}        \\\n'
                  '     -xLAB  {_xlab__:>12s}        \\\n'
                  '     -ndim  {_ndim__:>12s}        \\\n'
                  '     -out {_outputfile_} -verb -ov  \n'
                  '\n'
                  'nmrPipe   -in {_outputfile_} \\\n'
                  '| nmrPipe -fn POLY -time                                \\\n'
                  '| nmrPipe -fn SP -off 0.35 -end 0.95 -pow 2 -c 0.5      \\\n'
                  '| nmrPipe -fn ZF -auto -zf 2                            \\\n'
                  '| nmrPipe -fn FT                                        \\\n'
                  '| nmrPipe -fn PS -p0 {_p0x_:>6s} -p1 {_p1x_:>6s} -di -verb        \\\n'
                  '{Ext}| nmrPipe -fn EXT {ProtonExtraction:27s}           \\\n'
                  '| nmrPipe -fn POLY -auto                                \\\n'
                  '     -out {_processedfile_} -verb -ov\n').format(**parameters)
        return script
    ###################
    def script_for_regular_2D(self, parameters):
        """
        """
        script = ('#!/bin/csh\n'
                  '\n'
                  'var2pipe -in {_fidfile_} \\\n'
                  '     -xN    {_xN____:>12s}     -yN    {_yN____:>12s}    \\\n'
                  '     -xT    {_xT____:>12s}     -yT    {_yT____:>12s}    \\\n'
                  '     -xMODE {_xMODE_:>12s}     -yMODE {_yMODE_:>12s}    \\\n'
                  '     -xSW   {_sw____:>12s}     -ySW   {_sw2___:>12s}    \\\n'
                  '     -xOBS  {_xobs__:>12s}     -yOBS  {_yobs__:>12s}    \\\n'
                  '     -xCAR  {_XCAR__:>12s}     -yCAR  {_yCAR__:>12s}    \\\n'
                  '     -xLAB  {_xlab__:>12s}     -yLAB  {_ylab__:>12s}    \\\n'
                  '     -ndim  {_ndim__:>12s}     -aq2D  {_aq2D__:>12s}    \\\n'
                  '     -out {_outputfile_} -verb -ov                        \n'
                  '\n'
                  'nmrPipe   -in {_outputfile_} \\\n'
                  '| nmrPipe -fn POLY -time                                \\\n'
                  '| nmrPipe -fn SP -off 0.35 -end 0.95 -pow 2 -c 0.5      \\\n'
                  '| nmrPipe -fn ZF -auto -zf 2                            \\\n'
                  '| nmrPipe -fn FT                                        \\\n'
                  '| nmrPipe -fn PS -p0 {_p0x_:>6s} -p1 {_p1x_:>6s} -di -verb        \\\n'
                  '{Ext}| nmrPipe -fn EXT {ProtonExtraction:27s}           \\\n'
                  '| nmrPipe -fn TP                                        \\\n'
                  '{LP}| nmrPipe -fn LP -fb                                    \\\n'
                  '| nmrPipe -fn SP -off 0.35 -end 1.0 -pow 2 -c {_c_:5s}     \\\n'
                  '| nmrPipe -fn ZF -auto -zf 2                            \\\n'
                  '| nmrPipe -fn FT                                        \\\n'
                  '| nmrPipe -fn PS -p0 {_p0y_:>6s} -p1 {_p1y_:>6s} -di -verb        \\\n'
                  '{Rev}| nmrPipe -fn REV -sw                                   \\\n'
                  '| nmrPipe -fn TP                                        \\\n'
                  '| nmrPipe -fn POLY -auto                                \\\n'
                  '     -out {_processedfile_} -verb -ov\n').format(**parameters)
        return script
    ###################
    def script_for_semi_3D(self, parameters):
        """
        """
        script = ('#!/bin/csh\n'
                  '\n'
                  'var2pipe -in {_fidfile_} \\\n'
                  '     -xN    {_xN____:>12s}     -yN    {_eN____:>12s}    -zN    {_yN____:>12s}    \\\n'
                  '     -xT    {_xT____:>12s}     -yT    {_eT____:>12s}    -zT    {_yT____:>12s}    \\\n'
                  '     -xMODE {_xMODE_:>12s}     -yMODE {_eMODE_:>12s}    -zMODE {_yMODE_:>12s}    \\\n'
                  '     -xSW   {_sw____:>12s}     -ySW   {_esw2__:>12s}    -zSW   {_sw2___:>12s}    \\\n'
                  '     -xOBS  {_xobs__:>12s}     -yOBS  {_eobs__:>12s}    -zOBS  {_yobs__:>12s}    \\\n'
                  '     -xCAR  {_XCAR__:>12s}     -yCAR  {_eCAR__:>12s}    -zCAR  {_yCAR__:>12s}    \\\n'
                  '     -xLAB  {_xlab__:>12s}     -yLAB  {_elab__:>12s}    -zLAB  {_ylab__:>12s}    \\\n'
                  '     -ndim  {_ndim__:>12s}     -aq2D  {_aq2D__:>12s}                           \\\n'
                  '     -out {_outputfile_} -verb -ov                        \n'
                  '\n'
                  'xyz2pipe  -in {_outputfile_} \\\n'
                  '| nmrPipe -fn POLY -time                                \\\n'
                  '| nmrPipe -fn SP -off 0.35 -end 0.95 -pow 2 -c 0.5      \\\n'
                  '| nmrPipe -fn ZF -auto -zf 2                            \\\n'
                  '| nmrPipe -fn FT                                        \\\n'
                  '| nmrPipe -fn PS -p0 {_p0x_:>6s} -p1 {_p1x_:>6s} -di -verb        \\\n'
                  '{Ext}| nmrPipe -fn EXT {ProtonExtraction:27s}           \\\n'
                  '| nmrPipe -fn TP                                        \\\n'
                  '| nmrPipe -fn ZTP                                       \\\n'
                  '{LP}| nmrPipe -fn LP -fb                                    \\\n'
                  '| nmrPipe -fn SP -off 0.35 -end 1.0 -pow 2 -c {_c_:5s}     \\\n'
                  '| nmrPipe -fn ZF -auto -zf 2                            \\\n'
                  '| nmrPipe -fn FT                                        \\\n'
                  '| nmrPipe -fn PS -p0 {_p0y_:>6s} -p1 {_p1y_:>6s} -di -verb        \\\n'
                  '{Rev}| nmrPipe -fn REV -sw                                   \\\n'
                  '| nmrPipe -fn TP                                        \\\n'
                  '| nmrPipe -fn POLY -auto -verb                          \\\n'
                  '| pipe2xyz -x -out {_processedfile_} -verb -ov\n').format(**parameters)
        return script
    ###################
    def CreateConverFile(self,
                         userphase,
                         SecondDimension = None,
                         Fastprocess = None,
                         Extract = None,
                         Trosy_experiment = None):
        """
        It creates the convert file based
        Parameters:
        ===========
            * userphase = Phase 0 for proton phase correction
            * SecondDimension = 'N'itrogen or 'C'arbon for the proper second dim
            * Fastprocess = No linear prediction
            * Extract = extract proton dimension, like: [6.0,10.0] in ppm
            * Trosy_experiment = Trosy or non trosy
        Returns:
        ========
            * File created = self.ConvertFileName
        """
        # Pulse sequence
        print self.Info('seqfil')[0]+' pulse sequence was used'
        #
        if not SecondDimension or 'C' in SecondDimension:
            # Carbon
            SecDim = ['sw1','dfrq','dn']
        else:
            # Nitrogen
            SecDim = ['sw1','dfrq2','dn2']
        #
        if not Fastprocess:
            Fastprocess = False
        else:
            Fastprocess = True
        #
        parameters = {}
        parameters['_fidfile_'] = self.Path + self.FidFileName
        #
        xcar = self.Get_H2O_chemical_shift(float(self.temp))
        yobs = float(self.Info(SecDim[1])[0])
        phase = self.Get_Previous_Phase_values()

        ph0 = phase[0]
        try:
            phase[0] += float(userphase)
            phase[0] = phase[0] % 360.0
            print '----------------------------------------------------------------'
            print 'The value for p0 phase correction is {0:5.1f} + ('.format(ph0) + userphase + ') = {0:5.1f}'.format(phase[0])
            print '----------------------------------------------------------------'
        except ValueError:
            pass
        #
        # First dimension
        parameters['_xN____'] = self.Info('np')[0]
        parameters['_xT____'] = str(int(self.Info('np')[0]) / 2)
        parameters['_xMODE_'] = 'Complex'
        parameters['_sw____'] = '{0:10.4f}'.format(float(self.Info('sw')[0]))
        parameters['_xobs__'] = '{0:10.6f}'.format(float(self.Info('sfrq')[0]))
        parameters['_XCAR__'] = '{0:9.8f}'.format(xcar)
        parameters['_xlab__'] = self.Info('tn')[0]
        #
        # Second dimension
        parameters['_yN____'] = str(int(self.Info('ni')[0]) * 2)
        parameters['_yT____'] = self.Info('ni')[0]
        parameters['_yMODE_'] = 'Rance-Kay' #'Rance-Kay' or 'Complex'
        parameters['_sw2___'] = '{0:10.4f}'.format(float(self.Info(SecDim[0])[0]))
        parameters['_yobs__'] = '{0:10.6f}'.format(yobs)
        parameters['_yCAR__'] = '{0:9.6f}'.format(self.Get_carrier_in_PPM(xcar, float(self.Info('sfrq')[0]), yobs, SecondDimension))
        parameters['_ylab__'] = self.Info(SecDim[2])[0]
        #
        parameters['_ndim__'] = '2'
        parameters['_aq2D__'] = 'States'
        parameters['_outputfile_'] = self.Path + self.Get_Current_Dir() + '.fid'
        #
        parameters['_zf_'] = '-auto -zf 2'
        # Processing parameters
        if Fastprocess:
            parameters['LP'] = '#'
        else:
            parameters['LP'] = ''
        #
        parameters['_p0x_'] = '{0:5.1f}'.format(phase[0])
        parameters['_p1x_'] = '{0:5.1f}'.format(phase[1])
        parameters['_p0y_'] = '{0:5.1f}'.format(phase[2])
        parameters['_p1y_'] = '{0:5.1f}'.format(phase[3])
        #
        if not Extract:
            if SecondDimension == 'N' and self._2D:
                parameters['Ext'] = ''
            else:
                parameters['Ext'] = '#'
            parameters['ProtonExtraction'] = '-left -sw'
        else:
            parameters['Ext'] = ''
            parameters['ProtonExtraction'] = '-x1 {0:3.1f}ppm -xn {1:3.1f}ppm -sw -round 16'.format(max(Extract),min(Extract))
        #
        if self._2D and self.Info('f1180') and 'y' in self.Info('f1180')[0]:
            parameters['_c_']    = '{0:3.1f}'.format(1.0)
            parameters['_p0y_']  = '{0:5.1f}'.format(-90.0)
            parameters['_p1y_']  = '{0:5.1f}'.format(180.0)
        else:
            parameters['_c_']    = '{0:3.1f}'.format(0.5)
            parameters['_p0y_']  = '{0:5.1f}'.format(0.0)
            parameters['_p1y_']  = '{0:5.1f}'.format(0.0)
        #
        if Trosy_experiment:
            parameters['Rev'] = ''
        else:
            parameters['Rev'] = '#'
        #
        parameters['_processedfile_'] = self.Path + self.Get_Current_Dir() + '.dat'
        #
        #
        self.multiple_file = False # True if there are multiple files in the fid file
        #################
        if not self._2D:
            # 1D experiment
            parameters['_ndim__'] = '1'
            ################
            script = self.script_for_regular_1D(parameters)
            ################
        else:
            # 2D experiment
            self.onefile = 'single_hsqc'
            for element in self.Info('array')[0][1:-1].split(','):
                if not 'phase' in element:
                    self.onefile = element
            if self.onefile == 'single_hsqc':
                ################
                script = self.script_for_regular_2D(parameters)
                ################
            else:
                print 'NOTE: More than one measurements found in the fid file!'
                self.multiple_file = True # True if there are multiple files in the fid file
                parameters['_eN____'] = str(len(self.Info(self.onefile)[0].split()))
                parameters['_eT____'] = str(len(self.Info(self.onefile)[0].split()))
                parameters['_eMODE_'] = 'Real'
                parameters['_esw2__'] = '10.000'
                parameters['_eobs__'] = '10.000'
                parameters['_eCAR__'] = '5.000'
                parameters['_ndim__'] = '3'
                parameters['_elab__'] = '"Trel"'
                parameters['_outputfile_'] = self.Path + self.__temporary_folder + '/' + self.Get_Current_Dir() + '_%03d.fid'
                parameters['_processedfile_'] = self.Path + self.Get_Current_Dir() + '_%01d.dat'
                os.system('mkdir '+self.__temporary_folder)
                ################
                script = self.script_for_semi_3D(parameters)
                ################
        # Write out the convert file
        convertfile = open(self.ConvertFileName,'w')
        convertfile.write(script)
        convertfile.close()
        #
        print 'The script utilized to process the data:'
        print script
        return None
    ###################
    def RunConvertFile(self,
                       open_nmrDraw = True,
                       open_sparky = False,
                       nocleanup = False):
        """
        Run the convert script, show it in NmrDraw or Sparky and erase all
        files if needed
        Parameters:
        ===========
            * open_nmrDraw =
            * open_sparky =
            * nocleanup =
        """
        #
        os.system('chmod 755 ' + self.ConvertFileName)
        os.system('./' + self.ConvertFileName)
        #
        result_file = self.Path + self.Get_Current_Dir()
        if self.multiple_file:
            # If there are multiple file, the header conteins 3 dim => set to 2
            for i in range(len(self.Info(self.onefile)[0].split())):
                os.system('sethdr '+ result_file + '_' + str(i + 1) + '.dat -ndim 2')
                os.system('pipe2ucsf '+ result_file + '_' + str(i + 1) + '.dat '
                                      + result_file + '_' + str(i + 1) + '.ucsf')
        else:
            os.system('pipe2ucsf '+ result_file + '.dat '
                                  + result_file + '.ucsf')

        if open_nmrDraw:
            os.system('nmrDraw '+ result_file+ '.dat')

        if open_sparky:
            pass
            #os.system('sparky ' + result_file + '*.ucsf')

        if not nocleanup:
            if self.__temporary_folder in os.listdir('.'):
                os.system('rm -rf ' + self.__temporary_folder)
            os.system('rm -f ' + result_file + '*.fid')
            os.system('rm -f ' + result_file + '*.dat')
        #
        return None
    def ByeBye(self):
        """
        Write out the temperature to be sure the the carrier is at the right
        place!
        """
        warning = ''.join(('-'*60, '\n', '>'*10,
                           ' PLEASE MAKE SURE THAT YOU MEASURED @ ',
                           str(self.temp), ' C! ',
                           '<'*10, '\n', '-'*60, '\n'
                         ))
        sys.__stderr__.write(warning)
        return None
################################################################################



arguments = sys.argv

HC = Convert_HSQC(arguments)

