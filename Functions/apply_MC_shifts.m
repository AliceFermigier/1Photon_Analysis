%% Apply MC shifts on tif files
% written by Alice Fermigier @Herry lab 2025
% to be used for multichannel miniscope recordings

%#function isx

function [Movie_all, Max_all] = apply_MC_shifts(All_nam, All_nam_second,...
    T_DS_factor, Resize_factor, psf, Data_Format_in,...
    Chunk_size, Recording_Speed)

    % Correlation map / movie
    Max_all = cell(size(All_nam, 1), 1);
    Movie_all = cell(size(All_nam, 1), 1);

    for p = 1:size(All_nam, 1)
        % Identify Sessions of the same Animal
        reference_folder = All_nam{p}
        % Retrieve corresponding secondary channel data
        subfolder_extracted = %same as All-nam but G replaced by R
        
        % Correlation map
        Max = cell(size(subfolder_extracted, 2), 1);
        Movies = cell(size(subfolder_extracted, 2), 1);
        
        for pp = 1:size(reference_folder, 2)
            
            % Find all tif files in the folder
            current_reference_folder = [reference_folder{pp} '\'];
            Current_folder = %same as current_reference_folder but with R instead of G
            disp('Current_folder:')
            disp(Current_folder)
             
            % Get Folder Name
            [~,name_R,~] = fileparts(subfolder_extracted{pp});
            name = ['S' name_R];

            if ~exist([Current_folder 'processed_data'])
                mkdir([Current_folder 'processed_data'])
            else
                disp('File was already processed - Results will be overwritten !')
            end
            
            searchString = fullfile(Current_folder, ['*.' Data_Format_in]);
            Filenames = dir(searchString);
            disp('Filenames:')
            disp(Filenames)
            
            if strcmpi(Data_Format_in, 'tif') || strcmpi(Data_Format_in, 'tiff')
                [~, idx] = sort_nat({Filenames.date}); % Order the indices according to recording time
                
                num = size(Filenames, 1);
                disp(Filenames)
                disp(num)
                A = cell(num, 1);
                t = Tiff([Current_folder Filenames(idx(1)).name]);
                Size = size(read(t));
                
                % Resize data as indicated in the main script. The resize
                % factormust match the one used for reference channel
                
                try
                    parfor z = 1:num
                        A{z} = imresize(loadtiff([Current_folder Filenames(idx(z)).name]), ...
                            [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                    end
                catch
                    for z = 1:num
                        Temp_tif = loadtiff([Current_folder Filenames(idx(z)).name]);
                        
                        try
                            A{z} = imresize(Temp_tif, [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                        catch
                            if iscell(Temp_tif) && numel(Temp_tif) == 3
                                Piece1 = Temp_tif{1};
                                Piece2 = Temp_tif{3};                            
                                Temp_tif{2} = mean(cat(3, Piece1(:,:,end), Piece2(:,:,1)), 3, 'native');
                                A{z} = imresize(cat(3, Temp_tif{:}), [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                            else
                               error('This case has not occured before, please contact julian.hinz@fmi.ch if you cant fix it.')        
                            end
                        end
                    end
                end
                
                data = cat(3 , A{:});
                clear A
                Size = size(data);

                chunks_orig_s = [1:Chunk_size:ceil(Size(3)/Chunk_size)*Chunk_size];
                chunks_orig_e = [chunks_orig_s(2:end) - 1 Size(3)];
                chunks_orig_s(2:end) = chunks_orig_s(2:end) ;

                chunks_MCs = ceil(chunks_orig_s/T_DS_factor);                 
                chunks_MCe = ceil(chunks_orig_e/T_DS_factor); 

                numFiles = numel(chunks_orig_s);

                % Check the class of the tif files to make the conversion
                % into gray scale accurate
                Class_M = class(data);

                switch Class_M
                    case 'uint16'
                        Max_Value = 2^16 - 1;
                    case 'uint8'
                        Max_Value = 2^8 - 1;
                    otherwise
                         error('The data format is not uint, please check the format and run this script again !')
                end

                savefast([Current_folder 'processed_data\dataset'], 'data');
                Data_handle = matfile([Current_folder 'processed_data\dataset']);
                clear data

                MC = zeros(Size(1), Size(2), numel(1:T_DS_factor:Size(3)), Class_M);
                savefast([Current_folder 'processed_data\MC'], 'MC');
                clear MC
                MC_f = matfile([Current_folder 'processed_data\MC'], 'Writable', true);
            else
                disp('Data format not handled')
            end
            
            % Load MC shifts from the reference channel
            Shift_collection = load([current_reference_folder 'processed_data\MC_Shifts.mat'], 'Shift_collection');
            
            frame_counter = 0;

            for ppp = 1:numel(Shift_collection)

                Shifts = Shift_collection{ppp};   % [Nframes × 2]
                Nf = size(Shifts,1);

                for k = 1:Nf
                    frame_counter = frame_counter + 1;

                    frame = data(:,:,frame_counter);

                    dx = Shifts(k,1);
                    dy = Shifts(k,2);

                    frame_MC = imtranslate(frame, [dx dy], ...
                        'linear', 'FillValues', 0);

                    MC_f.MC(:,:,frame_counter) = frame_MC;
                end
            end
        end        
    end
end

