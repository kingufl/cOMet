cOMet

An error correction tool for optical mapping data

Quick steps:

    Compile the package in the root folder: make
    Go to bin folder
    Run cOMet: ./cOMet Rmap_input_file_path_and_filename

For running in parallel:

    Compile the package in the root folder: make parallel
    Go to bin folder
    Run cOMet: ./cOMet <Rmap_input_file_path_and_filename> <Number_of_parallel_streams> <stream_index> 
    
    
    For example, if total number of Rmaps is 1,000,000 and total number of streams is 1000. Then stream_index=1 will error correct Rmaps 1     to 1000. stream_index=2 will error correct Rmaps 1001 to 2000 and so on.
