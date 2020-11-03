To analyze data for the AnDi competion:

1. Run the code "make_directories.m" to create all required folders
2. Download and unzip the competition datasets (three files named 
task1.txt, task2.txt and task3.txt) from the AnDi website. 
3. Extract the trajectories using the extract_traj_task{X}.m files, where
X corresponds to the task number.
4. Run the corresponding taskX_Yd_skel.m file where X=task and Y=dim. 
On running a skel file, the output for each trajectory is saved in the 
corresponding folder inside the folder named "Submit".
5. The file "format_submit_task{X}.m" can then be used to combine the 
results into the format required for the AnDi challenge. 

