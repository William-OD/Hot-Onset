import Scripts.Generating_HEK as gen_hek
import Scripts.BackG_Onset_Processing as processing
import os

hek_dir = "HEK_Data/GOES15_HEK_Data.csv" #default HEK directory and filename read/write
output_dir = "Onset_Data/Onset_Data.csv" #default onset data directory and filename read/write

def generate_flare_database():
    print("You have selected to generate a flare database from the Heliophysics Event Knowledgebase (HEK).")
    print("This program will generate a database from HEK, which will then be used in the second program.")
    print("")
    data_dir = str(input("Please enter the data directory of the flares you wish to analyse: "))
    if not os.path.exists(data_dir):
        print("The specified data directory does not exist. Please make sure the directory is correct.")
        return
    print("Please note that this may take a few minutes.")
    print("The program will now begin to generate the database from HEK...")
    gen_hek.process_directory(data_dir, hek_dir)
    print("The program has finished generating the database from HEK.")
    print("The database has been saved to: " + hek_dir)
    print("")

def background_onset():
    print("You have selected to novel background subtract the data, and calculate the hot onset.")
    print("This program will novel background subtract the data, and calculate the hot onset data.")
    print("")
    data_dir = str(input("Please enter the data directory of the flares you wish to analyse: "))
    if not os.path.exists(data_dir):
        print("The specified data directory does not exist. Please make sure the directory is correct.")
        return
    print("")
    print("Please note that this may take a few hours, depending on the number of flares being analysed.")
    print("The program will now begin to novel background subtract the data, and calculate the hot onset data...")
    print("")
    default = str(input("Would you like to use the default settings? (y/n)"))
    print("")
    if default == "n":
        onset_length = int(input("Please enter the onset length in seconds (default = 60): "))
        print("")
        sigma = int(input("Please enter the onset sigma tolerance value (default = 4): "))
        print("")
        processing.process_data(data_dir, hek_dir, output_dir, onset_length, sigma)
    else:
        print("Using default settings...")
        processing.process_data(data_dir, hek_dir, output_dir, onset_length = 60, sigma = 4)
    print("The program has calculated the hot onset data and saved to: " + output_dir)

def both_programs():
    print("You have selected to generate a flare database from the Heliophysics Event Knowledgebase (HEK) and then carry out the hot onset analysis, after a novel background subtraction.")
    print("This will generate a database from HEK containing flare information from the earliest to latest flare in the specified directory, and saved. It will then novel background subtract the data, and extract the hot onsets.")
    print("")
    data_dir = str(input("Please enter the data directory of the flares you wish to analyse: "))
    if not os.path.exists(data_dir):
        print("The specified data directory does not exist. Please make sure the directory is correct.")
        return
    print("The first program will now begin to generate the database from HEK...")
    print("Please note that this may take a few minutes.")
    gen_hek.process_directory(data_dir, hek_dir)
    print("The program has finished generating the database from HEK.")
    print("The database has been saved to: " + hek_dir)
    print("")
    print("The second program will now begin to novel background subtract the data, and calculate the hot onset data...")
    print("Please note that this may take a few hours, depending on the number of flares being analysed.")
    print("")
    default = str(input("Would you like to use the default settings? (y/n)"))
    print("")
    if default == "n":
        onset_length = int(input("Please enter the onset length in seconds (default = 60): "))
        print("")
        sigma = int(input("Please enter the onset sigma tolerance value (default = 4): "))
        print("")
        processing.process_data(data_dir, hek_dir, output_dir, onset_length, sigma)
    else:
        print("Using default settings...")
        processing.process_data(data_dir, hek_dir, output_dir, onset_length = 60, sigma = 4)
    print("The program has calculated the hot onset data and saved to: " + output_dir)



print("")
print("------------------------------------------------------------")
print("Welcome to the Solar Flare Hot Onset Analysis Tool!")
print("")
print("Before use, please make sure that you have followed the correct setup instructions on the GitHub .Readme file.")
print("This includes following the correct environment setup, as well as downloading the correct data files.")
print("------------------------------------------------------------")
print("")

while True:
    print("")
    print("Which program would you like to run?")
    print("")
    print("1. Generate flare database from the Heliophysics Event Knowledgebase (HEK).")
    print("2. Calculating the hot onset data, after a novel background subtraction.")
    print("3. Run both programs 1 and 2.")
    print("")
    print("Type 'exit' to exit the program.")
    print("")
    # try:
    program_number = input("Please enter the number of the program you wish to run: ")
    print("")
    print("------------------------------------------------------------")
    print("")
    if program_number == "1":
        generate_flare_database()
        print("")
        print("------------------------------------------------------------")
    elif program_number == "2":
        background_onset()
    elif program_number == "3":
        print("You have selected to run both programs 1 and 2.")
        generate_flare_database()
        background_onset()
        print("Program 1 and 2 have been completed.")

    elif program_number == "exit":
        print("You have selected to exit the program.")
        print("Thank you for using the Solar Flare Hot Onset Analysis Tool! If you have any questions, please contact the author at: w.odonnell.1@research.gla.ac.uk")
        print("Exiting program...")
        print("")
        break

# TODO
# 1. Add in the option to change filename of output HEK file and output onset file (but not the directory).
# 2. Check over general customisation of the program. Any values that are non-fixed and could be expermineted with
# 3. Update in readme file (especially the setup instructions for downlaoding data and setting up repo). Outlining customisation options.
# 4. Add licences
# 5. Turn the GUI into a script.