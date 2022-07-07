classdef VevoRawReader < handle
    %A class to read and organize raw Vevo LAZR-X files

    properties
        Property1
    end
    
%     Define static methods
    methods (Static)
        function file_attributes = read_vevo_xml(filename)
            % reads the xml file into a cell array containing two columns: name and
            % value
            % input argument: filename (string)

            % read the xml header file
            try
                xDoc = xmlread(filename);
            catch
               error('Failed to read XML file %s.',filename);
            end
            % read elements by their tagname 'parameter'
            AllParameters = xDoc.getElementsByTagName('parameter');
            % loop through all the elements in the xml object and read
            % their name and value 
            file_attribute_name = {};
            file_attribute_value = {};
            for k = 0:AllParameters.getLength-1
                node = AllParameters.item(k);
                file_attribute_name{k+1} = char(node.getAttribute('name'));
                file_attribute_value{k+1} = char(node.getAttribute('value'));
            end
            % transpose the cell arrays to column form
            file_attribute_name = file_attribute_name';
            file_attribute_value = file_attribute_value';
            % concatenate the name and value arrays into a single cell
            % array
            file_attributes = [file_attribute_name , file_attribute_value];
        
        end

        function attribute_value = get_vevo_attribute(file_attributes,attribute_name)
            % reads the attribute value corresponding to a given attribute
            % name
            % input args: 
            % file_attributes (string): the output of the read_vevo_xml method
            % attribute_name (string): attribute_name of interest
            try
                % read the second column of file_attributes cell array (values)
                file_attribute_value = file_attributes(:,2);
                % find the attribute_name of interest and the corresponding
                % value
                attribute_value = file_attribute_value(string(file_attributes)==attribute_name);
                % convert to string
                attribute_value = char(attribute_value);
            catch
               error('Failed to find %s in header file.',attribute_name);
            end
        end
        
        function padded_image = padder(Rawdata,TopPad, BottomPad, LeftPad, RightPad)
            %PADDER Pads the Rawdata image on top, bottom, left and right
            % Rawdata: the 2D or 3D image of interest
            % TopPad, BottomPad, LeftPad, RightPad (numerical): pad values
            % in respective directions

            % first, pad the image on top and left directions
            RawDataLeftTopPadded = padarray(Rawdata,[TopPad, LeftPad], 'pre');
            % then, pad the image on bottom and right directions
            padded_image = padarray(RawDataLeftTopPadded,[BottomPad, RightPad], 'post');
        end
        
        function foldernamesFinal = read_Bruker_folders(dirName)
        
            list_of_folders = dir(dirName); %specify directory
            list_of_folders = list_of_folders([list_of_folders(:).isdir]); 
            foldernamesFinal = {list_of_folders.name}; %reading the file names
            foldernamesFinal = foldernamesFinal'; %transposing
            foldernamesFinal = foldernamesFinal(3:end); %the first two elements are . and .. (not wanted)
            foldernamesFinal = string(foldernamesFinal); %Convert to string array
            % [~,ind] = sort(str2double(foldernamesFinal));
            % foldernamesFinal = foldernamesFinal(ind);
            % In case we have some folders in the raw bruker files which do not have numbers as file name
            % foldernamesFinal(isnan(sort(str2double(foldernamesFinal)))) = [];
        end
    end



    methods
%         function obj = untitled4(inputArg1,inputArg2)
%             %UNTITLED4 Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end

        function scan_info = scan_information_table(obj, directory_of_interest)
            % creates a MATLAB table with information on scans in a given
            % directory
            % Input Args
            % directory_of_interest (string): location of the raw files
            % Output: MATLAB table with information on scans

            % find all files ending in .xml (header files)
            list_of_xmls = dir(fullfile(directory_of_interest, '*.raw.xml'));
            % extract the name field
            list_of_xmls = {list_of_xmls.name};
            % create empty fields to append to later in code
            subject_ids = {};
            mode_names = {};

            % Note that step_size information is stored even for 2D scans,
            % It is basically how the motor is set up for that particular
            % scan regardless of whether the scan type is 3D or not
            step_size_3d = {};
            total_steps_3d = {};
            scan_datetimes = {};
            directory_loc = {};
            % Loop through the xml header files
            for i = 1: length(list_of_xmls)
                filename = list_of_xmls{i};
                % Read the xml header into a file_attribute cell array
                file_attribute = obj.read_vevo_xml(fullfile(directory_of_interest,filename));
                % read several attributes from the header file
                subject_ids{i} = obj.get_vevo_attribute(file_attribute, 'Series-Name');
                step_size_3d{i} = obj.get_vevo_attribute(file_attribute, '3D-Step-Size');
                total_steps_3d{i} = str2double(obj.get_vevo_attribute(file_attribute, '3D-Scan-Distance')) / str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
                time_of_scan = obj.get_vevo_attribute(file_attribute, 'Acquired-Time');
                date_of_scan = obj.get_vevo_attribute(file_attribute, 'Acquired-Date');
                scan_datetime = append(date_of_scan,' ',time_of_scan);
                % Create a datetime object to allow sorting based on date
                % and time of scan
                scan_datetimes{i} = datetime(scan_datetime,'InputFormat','MM/dd/yyyy h:mm:ss a');
                directory_loc{i} = directory_of_interest;
                
                % Since OxyHemo scans have 'pamode' as their Mode-Name, we
                % distinguish them separately by detecting the file stored
                % instead
                if exist(fullfile(directory_of_interest, append(erase(filename,".raw.xml"),'.raw.3d.paoxy')),'file') ...
                        || exist(fullfile(directory_of_interest, append(erase(filename,".raw.xml"),'.raw.paoxy')),'file')
                    mode_names{i} = 'PaOxy-Mode';
                else
                    mode_names{i} = obj.get_vevo_attribute(file_attribute, 'Mode-Name');
                end

            end
            % Transpose the cell arrays to n*1 format
            subject_ids = subject_ids';
            mode_names = mode_names';
            step_size_3d = step_size_3d';
            total_steps_3d = total_steps_3d';
            scan_datetimes = scan_datetimes';

            % Erase the ".raw.xml" extension and save the filename
            file_names = erase(list_of_xmls,".raw.xml");
            file_names = file_names';
            % Convert the cell arrays into a single table containing all th
            %e scan information for a particular directory
            scan_info = cell2table([directory_loc, file_names, scan_datetimes, subject_ids,mode_names,step_size_3d,total_steps_3d],"VariableNames",["directory_loc", "file_names", "scan_datetimes", "subject_ids" "mode_names", "step_size_3d", "total_steps_3d"]);
%             Sort based on scan time
            scan_info = sortrows(scan_info,{'scan_datetimes'},{'ascend'});
        end
        

        function [Rawdata, WidthAxis, DepthAxis, ZAxis] = VsiOpenRawBmode8(obj,fnameBase)

            % VsiOpenRawBmode8.m
            % Copyright VisualSonics 1999-2012
            % A. Needles
            % Revision: 1.0 Oct 24 2012
            % A function to open RAW 8-bit B-Mode files from data export on the Vevo 2100
            % and read selected parameters
            
            % Updated by Ali Kamali Feb 2021

            % Rawdata is the raw bmode matrix 
            % WidthAxis and DepthAxis are from the bmode window
            % ZAxis depends on motor step size
            
            % fnameBase is the filename
            
            ModeName = '.bmode';
            
            % Set up file names
            file_attribute = obj.read_vevo_xml(append(fnameBase,'.raw.xml'));

            % Check to see if file is 3D
            Filename3D = [fnameBase '.raw' '.3d' ModeName]; %Create a 3D filename that Vevo outputs if it is a 3D file
            fstruct = dir(Filename3D); %Check to see if such file exists in the current directory
            if isempty(fstruct)
                fname = [fnameBase '.raw' ModeName]; % No .3d extension needed if the file is NOT 3D
                total_steps_3d = 1;
            else
                fname = Filename3D;
%                 total_steps_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Scan-Distance')) / str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
%                 total_steps_3d = int8(round(total_steps_3d));
            end
            
            fnameXml = [fnameBase '.raw' '.xml'];
            % Parse the XML parameter file - DO NOT CHANGE
            
            BmodeNumSamples = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Samples'));
            BmodeNumLines = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Lines'));
            BmodeDepthOffset = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Depth-Offset')); %mm
            BmodeDepth = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Depth')); %mm
            BmodeWidth = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Width')); %mm

%             param = VsiParseXml(fnameXml, '.pamode');
%             BmodeNumSamples = param.BmodeNumSamples;
%             BmodeNumLines = param.BmodeNumLines;
%             BmodeDepthOffset = param.BmodeDepthOffset; %mm
%             BmodeDepth =  param.BmodeDepth; %mm
%             BmodeWidth =  param.BmodeWidth; %mm
            step_size_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
            total_steps_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Scan-Distance')) / str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
            total_steps_3d = round(double(round(total_steps_3d)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            % This is to strip the header data in the files - DO NOT CHANGE
            size = 1; % - bytes
            file_header = 40; % 40bytes
            line_header = 0; % 4bytes
            frame_header = 56; % 56bytes  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fid = fopen(fname,'r');
            DepthAxis = [BmodeDepthOffset:(BmodeDepth-BmodeDepthOffset)/(BmodeNumSamples-1):BmodeDepth];
            WidthAxis = [0:BmodeWidth/(BmodeNumLines-1):BmodeWidth];
            ZAxis = [0:step_size_3d:(total_steps_3d-1)*step_size_3d];

            % Initialize data
            Rawdata = ones(BmodeNumSamples, BmodeNumLines,total_steps_3d,'int16');
            
            for iframe = 1:total_steps_3d
                header = file_header + frame_header*iframe + (size*BmodeNumSamples*BmodeNumLines + BmodeNumLines*line_header)*(iframe-1);
                for bmode_i=1:BmodeNumLines
                
                    fseek(fid, header + (size*BmodeNumSamples + line_header)*(bmode_i-1),-1);
                    fseek(fid, line_header, 'cof');
                    [Rawdata(:,bmode_i,iframe),~]=fread(fid, BmodeNumSamples, 'uchar');
                end    
            end
            
            fclose(fid);
        end

        function [RawdataOxy, RawdataHbT, MaskHbT, TopPad, BottomPad, LeftPad, RightPad, WidthAxis, DepthAxis, ZAxis] = open_raw_oxyhemo(obj,fnameBase, HbTThreshold)

            
            % VsiOpenRawOxyHemo.m
            
            % A function to open RAW OxyHemo-Mode files from data export on the Vevo
            % LAZR-X and read selected parameters
            
            
            % RawdataOxy is the raw Oxygenation image 
            % RawdataHbT is the raw total Hemoglobin image
            % MaskHbT is the HbT mask based on a given threshold 
            % e.g. HbT Threshold in the Vevolab software is set to "20" by default. 
            % What this means is that the lowest 20% of pixel values are thrown out, and the top 80% are used to generate a binary map.
            % The pad outputs specify the pad values necessary to pad the image into
            % bmode space (at pa resolution)
            % WidthAxis and DepthAxis are from the PAI window
            
            % fnameBase is the filename
            % iframe is the frame (or slice) number where data needs to be extracted
            % from
            % HbTThreshold
            ModeName = '.pamode';
            
            % Set up file names
            
            file_attribute = obj.read_vevo_xml(append(fnameBase,'.raw.xml'));
            % Check to see if file is 3D
            filename = [fnameBase '.raw' '.3d' ModeName]; %Create a 3D filename that Vevo outputs if it is a 3D file
            fstruct = dir(filename); %Check to see if such file exists in the current directory
            if isempty(fstruct)
                fname = [fnameBase '.raw' ModeName]; % No .3d extension needed if the file is NOT 3D
                total_steps_3d = 1;
            else
                fname = filename;
                total_steps_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Scan-Distance')) / str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
                total_steps_3d = int8(round(total_steps_3d));
            end
            
            PaNumSamples = str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Samples'));
            PaNumLines = str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Lines'));
            PaDepthOffset = str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Depth-Offset')); %mm
            PaDepth = str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Depth')); %mm
            PaWidth = str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Width')); %mm
            PaCentre =  - str2double(obj.get_vevo_attribute(file_attribute, 'Pa-Mode/Centre')); %mm
            BmodeNumSamples = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Samples'));
            BmodeNumLines = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Lines'));
            BmodeDepthOffset = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Depth-Offset')); %mm
            BmodeDepth = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Depth')); %mm
            BmodeWidth = str2double(obj.get_vevo_attribute(file_attribute, 'B-Mode/Width')); %mm

%             fnameXml = [fnameBase '.raw' '.xml'];
%             % Parse the XML parameter file - DO NOT CHANGE
%             param = VsiParseXml(fnameXml, '.paoxy');
%             PaNumSamples = param.PaNumSamples;
%             PaNumLines = param.PaNumLines;
%             PaDepthOffset = param.PaDepthOffset; %mm
%             PaDepth = param.PaDepth; %mm
%             PaWidth = param.PaWidth; %mm
%             PaCentre =  - param.PaCentre; %mm
%             BmodeNumSamples = param.BmodeNumSamples; %mm
%             BmodeNumLines = param.BmodeNumLines; %mm
%             BmodeDepthOffset = param.BmodeDepthOffset; %mm
%             BmodeDepth =  param.BmodeDepth; %mm
%             BmodeWidth =  param.BmodeWidth; %mm
            step_size_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
            total_steps_3d = str2double(obj.get_vevo_attribute(file_attribute, '3D-Scan-Distance')) / str2double(obj.get_vevo_attribute(file_attribute, '3D-Step-Size'));
            total_steps_3d = round(double(round(total_steps_3d)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % This is to strip the header data in the files - DO NOT CHANGE
            sizebyte = 2; % unsigned int bytes (The OxyHemo data is stored as unsigned integer which is 2 bytes.
            file_header = 40; % 40 bytes
            line_header = 0; % 0bytes no header
            frame_header = 56; % bytes  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Defining parameters to help with positioninng the color image inside the bmode
            % image. The formulas are directly copied from the C file.
            OffsetLeft = (PaCentre - PaWidth / 2) - (-BmodeWidth / 2);
            OffsetTop = PaDepthOffset - BmodeDepthOffset;
            
            
            % PA mode depth and width
            DepthAxis = [PaDepthOffset:(PaDepth-PaDepthOffset)/(PaNumSamples-1):PaDepth];
            WidthAxis = [OffsetLeft:(PaWidth-OffsetLeft)/(PaNumLines-1):PaWidth+OffsetLeft];
            ZAxis = [0:step_size_3d:(total_steps_3d-1)*step_size_3d];
            % WidthAxis = [0:PaWidth/(PaNumLines-1):PaWidth];
            % WidthAxis = WidthAxis - PaCentre/2;
            
            % Calculating the number of rows and columns to pad when putting in bmode
            % space. Rounding is because the division alomost never returns an integer.
            % The C code also truncates the values
            TopPad = round(OffsetTop/((PaDepth-PaDepthOffset)/(PaNumSamples-1)));
            LeftPad = round(OffsetLeft/(PaWidth/(PaNumLines-1)));
            BottomPad = round((BmodeDepth - max(DepthAxis))/((PaDepth-PaDepthOffset)/(PaNumSamples-1)));
            RightPad = round((BmodeWidth - max(WidthAxis))/(PaWidth/(PaNumLines-1)));
            
            
            
            
            % Initialize data
            % The OxyHemo data for oxygenation and HbT are stored subsequently per frame, hence two channels
            RawdataOxy = zeros(PaNumSamples, PaNumLines,total_steps_3d,'single');
            
            % First Channel (750nm data)
            fid = fopen(fname,'r');
            for iframe = 1:total_steps_3d
                header = file_header + frame_header*iframe + (2*sizebyte*PaNumSamples*PaNumLines + PaNumLines*line_header)*(iframe-1);
                for i=1:PaNumLines
                    fseek(fid, header + (sizebyte*PaNumSamples + line_header)*(i-1),-1);
                    fseek(fid, line_header, 'cof');
                    [RawdataOxy(:,i,iframe),~]=fread(fid, PaNumSamples, 'uint16'); %the PA data is store as unsigned int
                end
            end
            
            fclose(fid);
            
            % datatype range is between 0 : 256*256-1 (65535)
            RawdataOxy = 100*RawdataOxy/65535;
            
            
            % Second Channel (850nm data)
            % Initialize data
            RawdataHbT = zeros(PaNumSamples, PaNumLines,total_steps_3d,'uint16');
            
            fid = fopen(fname,'r');
            % The data from the second channel is stored right after the first channel
            % per frame. So skipping the binary file by PaNumSamples*PaNumLines
            % multiplied by the size (unsigned int = 2 bytes) takes us to the 850nm data
            % channel per frame.
            for iframe = 1:total_steps_3d
                header = sizebyte*PaNumSamples*PaNumLines + file_header + frame_header*iframe + (2*sizebyte*PaNumSamples*PaNumLines + PaNumLines*line_header)*(iframe-1);
                for i=1:PaNumLines
                    fseek(fid, header + (sizebyte*PaNumSamples + line_header)*(i-1),-1);
                    fseek(fid, line_header, 'cof');
                    [RawdataHbT(:,i,iframe),~]=fread(fid, PaNumSamples, 'uint16');
                end
            end
            fclose(fid);
            % Positioning the Hemo data in its bmode space
            % RawDataLeftTopPadded = padarray(RawdataHbT,[TopPad, LeftPad], 'pre');
            % RawdataHbTBmodeSpace = padarray(RawDataLeftTopPadded,[BottomPad, RightPad], 'post');
            
            % creating a mask based on the HbT thresholding
            MaskHbT = zeros(size(RawdataHbT),'logical');
            for iframe = 1:total_steps_3d
                RawdataHbT_2D= squeeze(RawdataHbT(:,:,iframe));
                MaskHbT_2D = RawdataHbT_2D(:);
                [~,TF] = rmoutliers(single(RawdataHbT_2D(:)),'percentiles',[HbTThreshold 100]);
                MaskHbT_2D(TF) = 0;
                MaskHbT_2D = reshape(MaskHbT_2D,size(RawdataHbT,1),size(RawdataHbT,2));
                MaskHbT_2D = MaskHbT_2D>0;
                MaskHbT(:,:,iframe) = MaskHbT_2D;
            end
            

        end
        
        

        function vevo_niftii_write_oxyhemo(obj, fnameBase, scanmode, output_filename, file_to_save, HbTThreshold)

            [RawdataOxy, RawdataHbT, MaskHbT, TopPad, BottomPad, LeftPad, RightPad, WidthAxis, DepthAxis, ZAxis] = open_raw_oxyhemo(obj,fnameBase, HbTThreshold);

            switch file_to_save
                case 'RawdataOxy'
                    ImageData3D = RawdataOxy;
                    ImageData3D(MaskHbT==0)=0;
                case 'RawdataHbT'
                    ImageData3D = RawdataHbT;
                case 'MaskHbT'
                    ImageData3D = MaskHbT;
            end
            
            % Bring into accompanying Bmode-Space
            ImageData3D = padder(ImageData3D,TopPad, BottomPad, LeftPad, RightPad);

            ImageData3D = rot90(rot90(ImageData3D));
            ImageData3D = permute(ImageData3D, [2,3,1]);
            ImageData3D = flip(ImageData3D,2);
            OriginalPixelDimensions = [DepthAxis(2) - DepthAxis(1), WidthAxis(2) - WidthAxis(1), ZAxis(2) - ZAxis(1)];
            data_type = class(ImageData3D);
            vevo_niftii_header = obj.make_niftii_header(scanmode, data_type, OriginalPixelDimensions, size(ImageData3D), output_filename);
            niftiwrite(ImageData3D,[output_filename '_' scanmode], vevo_niftii_header, 'Compressed',true);
        end

        function vevo_niftii_write_bmodeOH(obj, fnameBase, scanmode, output_filename)

            [Rawdata, WidthAxis, DepthAxis, ZAxis] = obj.VsiOpenRawBmode8(fnameBase);
            

            
            % Bring into accompanying Bmode-Space
            ImageData3D = rot90(rot90(Rawdata));
            ImageData3D = permute(ImageData3D, [2,3,1]);
            ImageData3D = flip(ImageData3D,2);
            OriginalPixelDimensions = [DepthAxis(2) - DepthAxis(1), WidthAxis(2) - WidthAxis(1), ZAxis(2) - ZAxis(1)];
            data_type = class(ImageData3D);
            vevo_niftii_header = obj.make_niftii_header(scanmode, data_type, OriginalPixelDimensions, size(ImageData3D), output_filename);
            niftiwrite(ImageData3D,[output_filename '_' scanmode], vevo_niftii_header, 'Compressed',true);
        end
        function vevo_niftii_header = make_niftii_header(obj,scanmode, data_type, voxel_dimensions_array, imagesize, output_filename)
            % initialize the struct file for the .nii header
            

            vevo_niftii_header = struct;
            % The next two lines follow the same logic behind changing the orientaion
            % of the image previously explained
            OriginalPixelDimensions = voxel_dimensions_array;
            vevo_niftii_header.PixelDimensions = [OriginalPixelDimensions(2), OriginalPixelDimensions(3),OriginalPixelDimensions(1)];
            
            % we will
            % need fields specific for a .nii file. The rest of this section adds these
            % the the header struct file
            vevo_niftii_header.Filename = output_filename;
            vevo_niftii_header.SpaceUnits = 'Millimeter';
            vevo_niftii_header.ImageSize = imagesize;
            vevo_niftii_header.Description = scanmode;
            vevo_niftii_header.Datatype = data_type;
            vevo_niftii_header.Description = 'None';
            vevo_niftii_header.Version = 'NIfTI1';
            vevo_niftii_header.Qfactor = 1;
            vevo_niftii_header.TimeUnits = 'None';
            vevo_niftii_header.SliceCode = 'Unknown';
            vevo_niftii_header.TransformName = 'Sform';
            
            % We should also position the image the right way when being read in
            % ITK-Snap, MIPAV etc.: So we make an affine3d transform object and include the pixel
            % dimensions in there. The transform type is similar to a T2-weighted MR
            % image
            
            tform = affine3d([vevo_niftii_header.PixelDimensions(1) 0 0 0; ...
                0 vevo_niftii_header.PixelDimensions(2) 0 0; ...
                0 0 vevo_niftii_header.PixelDimensions(3) 0; ...
                vevo_niftii_header.PixelDimensions(1) vevo_niftii_header.PixelDimensions(2) vevo_niftii_header.PixelDimensions(3) 1]);
            
            vevo_niftii_header.Transform = tform;
            vevo_niftii_header.FrequencyDimension = 0;
            vevo_niftii_header.PhaseDimension = 0;
            vevo_niftii_header.SpatialDimension = 0;


        end


        




    end
end