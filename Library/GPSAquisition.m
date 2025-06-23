classdef GPSAquisition < matlab.System
    % GPSAquisition encapsulates GNSS signal acquisition and processing
    % based on Kai SoftGNSS Book.

    % Public, non-tunable properties
    % properties (Nontunable)
    properties (Nontunable)
         
        %% Raw Signal parameter ===============================
        
        % Is Input Complex?
        isInputComplex (1,1) logical = true;
        % Sampling Frequency [Hz]
        samplingFreq = 2.4e6;   
        % Intermediate Frequency [Hz]
        IF = 0e6;  

         %% Processing settings ====================================================

        % Number of channels to be used for signal processing
         numberOfChannels = 8; 
 
        %% Acquisition settings ===================================================
        
        % List of satellites to look for. Some satellites can be excluded to speed up acquisition
        %[PRN numbers]PRN32 was just launched (as of 2014 Apr 15) so it was not set active yet:accuracy indicated by decoded message:11
        acqSatelliteList = [1:29 31 32];    
        % Band around IF to search for satellite signal. Depends on max Doppler [kHz]
        acqSearchBand = 14; 
        % Threshold for the signal presence decision rule
         % acqThreshold = 2.5;
         acqThreshold = 3.0;
        
        % No. of code periods for coherent integration (multiple of 2)
        acquisitionCohCodePeriods = 2;  
        % No. of non-coherent summations
        acquisitionNonCohSums = 4;     

        %% Constants ==============================================================
        
        % Initial sign. travel time  [ms] 
        % startOffset = 68.802;   

        startOffset = 0;   

        %% Block Related ==============================================================

        % Output Sampling time of the Block?
        UseSampleTimeOutput (1,1) logical = true;
        % Output Toggle/done of the Block?
        UseToggleOutput (1,1) logical = true;
        % Sample Time (-1 for inherited), [suggested >=0.25 sec]
        SampleTime = 0.5;

        

        % var1 (1,1) logical = false      % check box
        % var2 (1,1) {mustBePositive, mustBeInteger} = 10    % input box
        % var3 (1,:) char {mustBeMember(var3, ['opt1','opt2','opt3'])} = 'opt2' % drop down / combo box   

    end

    properties (Constant, Access = private)
        validRange = [1, 32]; % Valid range for acqSatelliteList elements
    end

    % properties (Hidden)
    % 
    % 
    % end
    % Public, tunable properties
    % properties
    % 
    % end

    % Pre-computed constants or internal states
    properties (Access = private)
        previousPeaks;
        fileType;
        caCodesTable;
        NoSamps;
        toggle;
        
        codeFreqBasis = 1.023e6;    % GPS stuff, constant
        codeLength = 1023; % GPS stuff, constant
        c = 299792458;  % Light speed, constant
        acquisitionTimeMS = 11;  % Dont change
        % acquisitionTimeMS = 30;



    end
  

    methods (Access = protected)

        function validatePropertiesImpl(obj)
            % Validate and initialize properties that depend on others

            % Validate isInputComplex
            validateattributes(obj.isInputComplex, {'logical'}, {'scalar', 'nonempty'});
            
            % Validate samplingFreq (positive scalar)
            validateattributes(obj.samplingFreq, {'numeric'}, {'positive', 'scalar', 'nonempty'});
            
            % Validate IF (non-negative scalar)
            validateattributes(obj.IF, {'numeric'}, {'nonnegative', 'scalar', 'nonempty'});
            
            % Validate numberOfChannels (positive integer scalar)
            validateattributes(obj.numberOfChannels, {'numeric'}, {'positive', 'integer', 'scalar', 'nonempty'});
            
            % Validate acqSatelliteList (vector of integers within range [1, 32])
            validateattributes(obj.acqSatelliteList, {'numeric'}, {'vector', 'integer', 'nonempty'});
            mustBeInRange(obj.acqSatelliteList, 1, 32);
            
            % Validate acqSearchBand (positive scalar)
            validateattributes(obj.acqSearchBand, {'numeric'}, {'positive', 'scalar', 'nonempty'});
            
            % Validate acqThreshold (positive scalar)
            validateattributes(obj.acqThreshold, {'numeric'}, {'positive', 'scalar', 'nonempty'});
            
            % Validate acquisitionCohCodePeriods (positive even integer)
            validateattributes(obj.acquisitionCohCodePeriods, {'numeric'}, {'positive', 'integer', 'scalar', 'nonempty'});
            assert(mod(obj.acquisitionCohCodePeriods, 2) == 0, 'acquisitionCohCodePeriods must be an even number.');
            
            % Validate acquisitionNonCohSums (positive integer scalar)
            validateattributes(obj.acquisitionNonCohSums, {'numeric'}, {'positive', 'integer', 'scalar', 'nonempty'});
            
            % Validate SampleTime (scalar, -1 for inherited or non-negative)
            validateattributes(obj.SampleTime, {'numeric'}, {'scalar', 'nonempty'});
            assert(obj.SampleTime == -1 || obj.SampleTime >= 0.25, 'SampleTime must be either -1 or >= 0.25.');


            % Initialize obj.NoSamps early
            obj.NoSamps = round(obj.samplingFreq * obj.acquisitionTimeMS * 1e-3);            
        end
        
        function setupImpl(obj, ~)
            % Perform one-time calculations, such as computing constants
            % and Initialization

            % Initialize fileType based on isInputComplex
            if obj.isInputComplex
                obj.fileType = 2;
            else
                obj.fileType = 1;
            end
            obj.previousPeaks = zeros(32, 1);
            obj.caCodesTable = obj.makeCaTable();
            obj.toggle = 0;
        end

        
        %%% Main Function %%%
        function [satPeaks, acqThreshold, varargout] = stepImpl(obj, samples)
        
            % Ensure the data is double before complex arithmetic, since
            % complex arithmetic is not supported for other types like integers.
            samples = double(samples);  
            
            % Ensure samples is a row vector
            samples = samples(:)'; % Acquisition assumes samples are provided as a row vector.
                  
            % Perform acquisition
            acqResults = obj.acquisition(samples, obj.acqThreshold);
            satPeaks = acqResults.peakMetric;
            satPeaks = satPeaks(:);  % Ensure return is a column vector
            acqThreshold = obj.acqThreshold;
        
            % If for some reason any element in satPeaks is NaN, utilize previous peaks
            if any(isnan(satPeaks))
                satPeaks = obj.previousPeaks;
            else
                obj.previousPeaks = satPeaks;
            end
        
            % Initialize varargout index
            varargoutIndex = 1;
            
            % Additional outputs based on object properties
            if obj.UseSampleTimeOutput
                sts = getSampleTime(obj);
                varargout{varargoutIndex} = sts.SampleTime;
                varargoutIndex = varargoutIndex + 1;
            end
        
            if obj.UseToggleOutput
                varargout{varargoutIndex} = obj.toggle;
                obj.toggle = double(~obj.toggle); % Toggle the value for the next step
            end
        end




        
   

        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% My functions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate CA Code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function CAcode = generateCAcode(~, PRN)
        % generateCAcode.m generates one of the 32 GPS satellite C/A codes.
        %
        % CAcode = generateCAcode(PRN)
        %
        %   Inputs:
        %       PRN         - PRN number of the sequence.
        %
        %   Outputs:
        %       CAcode      - a vector containing the desired C/A code sequence 
        %                   (chips).  
        
        %--------------------------------------------------------------------------

        %CVS record:
        %$Id: generateCAcode.m,v 1.1.2.5 2006/08/14 11:38:22 dpl Exp $
        
        %--- Make the code shift array. The shift depends on the PRN number -------
        % The g2s vector holds the appropriate shift of the g2 code to generate
        % the C/A code (ex. for SV#19 - use a G2 shift of g2s(19) = 471)
        g2s = [  5,   6,   7,   8,  17,  18, 139, 140, 141, 251, ...
               252, 254, 255, 256, 257, 258, 469, 470, 471, 472, ...
               473, 474, 509, 512, 513, 514, 515, 516, 859, 860, ...
               861, 862 ... end of shifts for GPS satellites 
               ... Shifts for the ground GPS transmitter are not included
               ... Shifts for EGNOS and WAAS satellites (true_PRN = PRN + 87)
                         145, 175,  52,  21, 237, 235, 886, 657, ...
               634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, ...
               386];
        
        %--- Pick right shift for the given PRN number ----------------------------
        g2shift = g2s(PRN);
        
        %--- Generate G1 code -----------------------------------------------------
        
        %--- Initialize g1 output to speed up the function ---
        g1 = zeros(1, 1023);
        %--- Load shift register ---
        reg = -1*ones(1, 10);
        
        %--- Generate all G1 signal chips based on the G1 feedback polynomial -----
        for i=1:1023
            g1(i)       = reg(10);
            saveBit     = reg(3)*reg(10);
            reg(2:10)   = reg(1:9);
            reg(1)      = saveBit;
        end
        
        %--- Generate G2 code -----------------------------------------------------
        
        %--- Initialize g2 output to speed up the function ---
        g2 = zeros(1, 1023);
        %--- Load shift register ---
        reg = -1*ones(1, 10);
        
        %--- Generate all G2 signal chips based on the G2 feedback polynomial -----
        for i=1:1023
            g2(i)       = reg(10);
            saveBit     = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10);
            reg(2:10)   = reg(1:9);
            reg(1)      = saveBit;
        end
        
        %--- Shift G2 code --------------------------------------------------------
        %The idea: g2 = concatenate[ g2_right_part, g2_left_part ];
        g2 = [g2(1023-g2shift+1 : 1023), g2(1 : 1023-g2shift)];
        
        %--- Form single sample C/A code by multiplying G1 and G2 -----------------
        CAcode = -(g1 .* g2);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make CA Tables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function caCodesTable = makeCaTable(obj)
        %Function generates CA codes for all 32 satellites based on the settings
        %provided in the structure "settings". The codes are digitized at the
        %sampling frequency specified in the settings structure.
        %One row in the "caCodesTable" is one C/A code. The row number is the PRN
        %number of the C/A code.
        %
        %caCodesTable = makeCaTable(settings)
        %
        %   Inputs:
        %       settings        - receiver settings
        %   Outputs:
        %       caCodesTable    - an array of arrays (matrix) containing C/A codes
        %                       for all satellite PRN-s
        
        %--------------------------------------------------------------------------

        
        %CVS record:
        %$Id: makeCaTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $
        
        %--- Find number of samples per spreading code ----------------------------
        samplesPerCode = round(obj.samplingFreq / ...
                                   (obj.codeFreqBasis / obj.codeLength));
        
        %--- Prepare the output matrix to speed up function -----------------------
        caCodesTable = zeros(32, samplesPerCode);
         
        %--- Find time constants --------------------------------------------------
        ts = 1/obj.samplingFreq;   % Sampling period in sec
        tc = 1/obj.codeFreqBasis;  % C/A chip period in sec
         
        %=== For all satellite PRN-s ...
        for PRN = 1:32
            %--- Generate CA code for given PRN -----------------------------------
            caCode = obj.generateCAcode(PRN);
         
            %=== Digitizing =======================================================
            
            %--- Make index array to read C/A code values -------------------------
            % The length of the index array depends on the sampling frequency -
            % number of samples per millisecond (because one C/A code period is one
            % millisecond).
            codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
            
            %--- Correct the last index (due to number rounding issues) -----------
            codeValueIndex(end) = 1023;
            
            %--- Make the digitized version of the C/A code -----------------------
            % The "upsampled" code is made by selecting values form the CA code
            % chip array (caCode) for the time instances of each sample.
            caCodesTable(PRN, :) = caCode(codeValueIndex);
            
        end % for PRN = 1:32
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Acquisition  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %-----------------------------------------------------------------------------------
        %-----------------------------------------------------------------------------------
        function acqResults = acquisition(obj, longSignal, acqThreshold)
        %Function performs cold start acquisition on the collected "data". It
        %searches for GPS signals of all satellites, which are listed in field
        %"acqSatelliteList" in the settings structure. Function saves code phase
        %and frequency of the detected signals in the "acqResults" structure.
        %
        %acqResults = acquisition(longSignal, settings)
        %
        %   Inputs:
        %       longSignal    - 11 ms of raw signal from the front-end 
        %       settings      - Receiver  Provides information about
        %                       sampling and intermediate frequencies and other
        %                       parameters including the list of the satellites to
        %                       be acquired.
        %   Outputs:
        %       acqResults    - Function saves code phases and frequencies of the 
        %                       detected signals in the "acqResults" structure. The
        %                       field "carrFreq" is set to 0 if the signal is not
        %                       detected for the given PRN number. 
         
        %--------------------------------------------------------------------------
        
        %CVS record:
        %$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $
        
        %% Initialization =========================================================
        
        % Find number of samples per spreading code
        samplesPerCode = round(obj.samplingFreq / ...
                                (obj.codeFreqBasis / obj.codeLength));
        
        % Create two 1msec vectors of data to correlate with and one with zero DC
        signal1 = longSignal(1 : samplesPerCode);
        signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);
        
        signal0DC = longSignal - mean(longSignal); 
        
        % Find sampling period
        ts = 1 / obj.samplingFreq;
        
        % Find phase points of the local carrier wave 
        phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;
        
        % Number of the frequency bins for the given acquisition band (500Hz steps)
        numberOfFrqBins = round(obj.acqSearchBand * 2) + 1;
        
        % % Generate all C/A codes and sample them according to the sampling freq.
        % caCodesTable = makeCaTable(settings);
        
        
        %--- Initialize arrays to speed up the code -------------------------------
        % Search results of all frequency bins and code shifts (for one satellite)
        results     = zeros(numberOfFrqBins, samplesPerCode);
        
        % Carrier frequencies of the frequency bins
        frqBins     = zeros(1, numberOfFrqBins);
        
        
        %--- Initialize acqResults ------------------------------------------------
        % Carrier frequencies of detected signals
        acqResults.carrFreq     = zeros(1, 32);
        % C/A code phases of detected signals
        acqResults.codePhase    = zeros(1, 32);
        % Correlation peak ratios of the detected signals
        acqResults.peakMetric   = zeros(1, 32);
        
        % fprintf('(');
        
        % Perform search for all listed PRN numbers ...
        for PRN = obj.acqSatelliteList
        
        %% Correlate signals ======================================================   
            %--- Perform DFT of C/A code ------------------------------------------
            caCodeFreqDom = conj(fft(obj.caCodesTable(PRN, :)));
        
            %--- Make the correlation for whole frequency band (for all freq. bins)
            for frqBinIndex = 1:numberOfFrqBins
        
                %--- Generate carrier wave frequency grid (0.5kHz step) -----------
                frqBins(frqBinIndex) = obj.IF - ...
                                       (obj.acqSearchBand/2) * 1000 + ...
                                       0.5e3 * (frqBinIndex - 1);
        
                %--- Generate local sine and cosine -------------------------------
                sigCarr = exp(1i*frqBins(frqBinIndex) * phasePoints);
                
                %--- "Remove carrier" from the signal -----------------------------
                sig1      = sigCarr .* signal1;
                sig2      = sigCarr .* signal2;

                %--- Convert the baseband signal to frequency domain --------------
                IQfreqDom1 = fft(sig1);
                IQfreqDom2 = fft(sig2);

                                % %--- "Remove carrier" from the signal -----------------------------
                % I1      = real(sigCarr .* signal1);
                % Q1      = imag(sigCarr .* signal1);
                % I2      = real(sigCarr .* signal2);
                % Q2      = imag(sigCarr .* signal2);
                % 
                % %--- Convert the baseband signal to frequency domain --------------
                % IQfreqDom1 = fft(I1 + 1i*Q1);
                % IQfreqDom2 = fft(I2 + 1i*Q2);
        
        
                %--- Multiplication in the frequency domain (correlation in time
                %domain)
                convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
                convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
        
                %--- Perform inverse DFT and store correlation results ------------
                acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
                acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
                
                %--- Check which msec had the greater power and save that, will
                %"blend" 1st and 2nd msec but will correct data bit issues
                if (max(acqRes1) > max(acqRes2))
                    results(frqBinIndex, :) = acqRes1;
                else
                    results(frqBinIndex, :) = acqRes2;
                end
                
            end % frqBinIndex = 1:numberOfFrqBins
        
        %% Look for correlation peaks in the results ==============================
            % Find the highest peak and compare it to the second highest peak
            % The second peak is chosen not closer than 1 chip to the highest peak
            
            %--- Find the correlation peak and the carrier frequency --------------
            [peakSize, frequencyBinIndex] = max(max(results, [], 2));
        
            %--- Find code phase of the same correlation peak ---------------------
            [peakSize, codePhase] = max(max(results));
        
            %--- Find 1 chip wide C/A code phase exclude range around the peak ----
            samplesPerCodeChip   = round(obj.samplingFreq / obj.codeFreqBasis);
            excludeRangeIndex1 = codePhase - samplesPerCodeChip;
            excludeRangeIndex2 = codePhase + samplesPerCodeChip;
        
            %--- Correct C/A code phase exclude range if the range includes array
            %boundaries
            if excludeRangeIndex1 < 2
                codePhaseRange = excludeRangeIndex2 : ...
                                 (samplesPerCode + excludeRangeIndex1);
                                 
            elseif excludeRangeIndex2 >= samplesPerCode
                codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                                 excludeRangeIndex1;
            else
                codePhaseRange = [1:excludeRangeIndex1, ...
                                  excludeRangeIndex2 : samplesPerCode];
            end
        
            %--- Find the second highest correlation peak in the same freq. bin ---
            % Ensure frequencyBinIndex is within valid range (added by me)
            if frequencyBinIndex < 1 || frequencyBinIndex > size(results, 1)
                disp('Frequency bin index out of range: ');
            else
                % Ensure codePhaseRange is within valid range
                codePhaseRange = codePhaseRange(codePhaseRange >= 1 & codePhaseRange <= size(results, 2));
                if isempty(codePhaseRange)
                    disp('Code phase range out of valid range.');
                else
                    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
                    %--- Store result -----------------------------------------------------
                    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
                    
                    % If the result is above threshold, then there is a signal ...
                    if (peakSize/secondPeakSize) > acqThreshold
                
                    %% Fine resolution frequency search =======================================
                        
                        %--- Indicate PRN number of the detected signal -------------------
                        % fprintf('%02d ', PRN);
                        
                        %--- Generate 10msec long C/A codes sequence for given PRN --------
                        caCode = obj.generateCAcode(PRN);
                        
                        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                                               (1/obj.codeFreqBasis));
                                           
                        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
                    
                        %--- Remove C/A code modulation from the original signal ----------
                        % (Using detected C/A code phase)
                        xCarrier = ...
                            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
                            .* longCaCode;
                        
                        %--- Compute the magnitude of the FFT, find maximum and the
                        %associated carrier frequency
                        
                        %--- Find the next highest power of two and increase by 8x --------
                        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
                        
                        %--- Compute the magnitude of the FFT, find maximum and the
                        %associated carrier frequency 
                        fftxc = abs(fft(xCarrier, fftNumPts)); 
                        
                        
                        uniqFftPts = ceil((fftNumPts + 1) / 2);
                        [fftMax, fftMaxIndex] = max(fftxc);
                        fftFreqBins = (0 : uniqFftPts-1) * obj.samplingFreq/fftNumPts;
                        if (fftMaxIndex > uniqFftPts) %and should validate using complex data
                            if (rem(fftNumPts,2)==0)  %even number of points, so DC and Fs/2 computed
                                fftFreqBinsRev=-fftFreqBins((uniqFftPts-1):-1:2);
                                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                                acqResults.carrFreq(PRN)  = -fftFreqBinsRev(fftMaxIndex);
                            else  %odd points so only DC is not included
                                fftFreqBinsRev=-fftFreqBins((uniqFftPts):-1:2);
                                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                                acqResults.carrFreq(PRN)  = fftFreqBinsRev(fftMaxIndex);
                            end
                        else
                            acqResults.carrFreq(PRN)  = (-1)^(obj.fileType-1)*fftFreqBins(fftMaxIndex);
                        end
                        
                        acqResults.codePhase(PRN) = codePhase;
                    
                    end   % if (peakSize/secondPeakSize) > acqThreshold
                    
                end
            end
        
        end    % for PRN = satelliteList
        end
        %=== Acquisition is over ==================================================
    end

%%%%%%%%%%%%%%%% Simulink Block Related Functions %%%%%%%%%%%%%%%% 
    methods (Access=protected)
        

        function resetImpl(obj)
            % Initialize / reset internal properties
            obj.toggle = 0;

        end

        function sts = getSampleTimeImpl(obj)
               if obj.SampleTime == -1
                       sts = createSampleTime(obj,'Type','Inherited');
               else
                       sts = createSampleTime(obj,'Type','Discrete',...
                         'SampleTime',obj.SampleTime);
               end
        end

        function validateInputsImpl(obj, samples)
            % Validate the size of the input vector samples
        
            % Validate samples:
            validateattributes(samples, {'numeric'}, {'vector', 'numel', obj.NoSamps}, '', 'samples');
        end

        function num = getNumInputsImpl(~)
            % Start with the base number of inputs (1 mandatory input)
            num = 1;  % Base input is "Samples"
        end

        function varargout = getInputNamesImpl(~)
            % Return the names of the inputs
            varargout = {'Samples'};  % Base input name
        end

        function num = getNumOutputsImpl(obj)
            % Start with the base number of outputs (2 mandatory outputs)
            num = 2;

            % Check if the optional output for sample time is enabled
            if obj.UseSampleTimeOutput
                num = num + 1;  % Add one more output
            end

            % Check if the optional output for toggle is enabled
            if obj.UseToggleOutput
                num = num + 1;  % Add one more output
            end

        end

        function varargout = getOutputNamesImpl(obj)
            % Return the names of the outputs
            varargout = {'Acqn Metric', 'Acqn Threshold'};  % Base output name

            if obj.UseSampleTimeOutput
                varargout{end+1} = 'Sample Time';
            end

            if obj.UseToggleOutput
                varargout{end+1} = 'Toggle/Done';
            end

        end
        




        % function idx = getSampleTimeOuputIdx(obj)
        %     % Get index of SampleTime output if enabled
        %     if obj.UseSampleTimeOutput
        %         idx = 3;
        %     else
        %         idx = [];
        %     end
        % end
        % 
        % 
        % function idx = getToggleOutputIdx(obj)
        %     % Get index of toggle output if enabled
        %     if obj.UseSampleTimeOutput && obj.UseToggleOutput
        %         idx = 4;
        %     elseif obj.UseToggleOutput
        %         idx = 3;
        %     else
        %         idx = [];
        %     end
        % end
        % An alternate approch
        % function idx = getOutputIdx(obj, outputType)
        %     if strcmp(outputType, 'SampleTime') && obj.UseSampleTimeOutput
        %         idx = 3;
        %     elseif strcmp(outputType, 'Toggle') && obj.UseToggleOutput
        %         idx = obj.UseSampleTimeOutput ? 4 : 3;
        %     else
        %         idx = [];
        %     end
        % end

      

        function varargout = getOutputSizeImpl(obj)
            peakSize = [32, 1];
            scalarSize = [1, 1];
        
            % Base outputs
            varargout{1} = peakSize;  % Peaks
            varargout{2} = scalarSize;  % Acquisition Threshold
        
            % Conditional outputs
            outputIndex = 3;
            if obj.UseSampleTimeOutput
                varargout{outputIndex} = scalarSize;  % Time Taken by loop / sampling time
                outputIndex = outputIndex + 1;
            end
            if obj.UseToggleOutput
                varargout{outputIndex} = scalarSize;  % Toggle
            end
        end



        function varargout = getOutputDataTypeImpl(obj)
            % Base outputs
            varargout{1} = "double";  % Peaks
            varargout{2} = "double";  % Acquisition Threshold
        
            % Conditional outputs
            outputIndex = 3;
            if obj.UseSampleTimeOutput
                varargout{outputIndex} = "double";  % Time Taken by loop / sampling time
                outputIndex = outputIndex + 1;
            end
            if obj.UseToggleOutput
                varargout{outputIndex} = "double";  % Toggle
            end
        end


        function varargout = isOutputComplexImpl(obj)
            % Base outputs
            varargout{1} = false;  % Peaks
            varargout{2} = false;  % Acquisition Threshold
        
            % Conditional outputs
            outputIndex = 3;
            if obj.UseSampleTimeOutput
                varargout{outputIndex} = false;  % Time Taken by loop / sampling time
                outputIndex = outputIndex + 1;
            end
            if obj.UseToggleOutput
                varargout{outputIndex} = false;  % Toggle
            end
        end


        function varargout = isOutputFixedSizeImpl(obj)
            % Base outputs
            varargout{1} = true;  % Peaks
            varargout{2} = true;  % Acquisition Threshold
        
            % Conditional outputs
            outputIndex = 3;
            if obj.UseSampleTimeOutput
                varargout{outputIndex} = true;  % Time Taken by loop / sampling time
                outputIndex = outputIndex + 1;
            end
            if obj.UseToggleOutput
                varargout{outputIndex} = true;  % Toggle
            end
        end


        function str = getIconImpl(~)
            str = sprintf('GPS\nAcquisition');
        end   
    end

    methods (Static, Access = protected)    
        
        function header = getHeaderImpl()
            % Define block header information
            header = matlab.system.display.Header(mfilename('class'), ...
                'Title', 'GPSAcquisition', ...
                'Text', 'GPSAquisition encapsulates GNSS signal acquisition and processing based on SoftGNSS Book.');
        end

        
        % function simMode = getSimulateUsingImpl
        %     % Return only allowed simulation mode in Simulink
        %     simMode = 'Code generation';
        % end
        
        function flag = showSimulateUsingImpl
            % Always show simulation mode combo box
            flag = true;
        end
    
    end


     methods
        function carrFreq = getDoppler(obj, samples)
            % getDoppler  Return per‐PRN Doppler estimates (Hz)
            %   samples - raw 11 ms IQ vector (row or column)
            %   carrFreq - 32×1 vector of estimated carrier freqs from acquisition()

            % Ensure double row vector
            samples = double(samples(:).');
            % Run the internal acquisition routine to get full struct
            acqResults = obj.acquisition(samples, obj.acqThreshold);
            % Return the Doppler for all 32 PRNs
            carrFreq = acqResults.carrFreq(:);
        end
    end

end
