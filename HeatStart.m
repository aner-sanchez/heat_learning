[~,current_dir,~] = fileparts(pwd);
if current_dir == "heat_learning"
    if exist('cvx', 'dir') == 7
        cd cvx\
        if exist('cvx_setup','file') == 2
            cvx_setup
            cd ..
            cd aux_functions\
            addpath(pwd)
            cd ..
            clc
            fprintf("All ready!\n" + ...
                "____________________________________\n" + ...
                "Type DemoHeat in the terminal for a demonstration of what I did\n\n")
        else
            fprintf("You have a directory called cvx but not a file called cvx_setup" + ...
                "it is strange")
        end
    else
        fprintf("You have to go here: http://cvxr.com/cvx/download/\n" + ...
            "to download cvx, it is a package for convex optimization" + ...
            "I use it for a function which I can't get myself");
        fprintf("..................................\n");
        fprintf("You need to download it and unzip it in heat_learning/cvx")
    end
else
    disp("go to directory ~\heat_learning and type it again");
end