% Iñigo Basterretxea Jacob
% 01246662
% This file prompts the user to enter the  
clear all;
close all;
clc;

Lx = input('Enter width of structure\n');
Ly = input('Enter height of structure\n');

nx = input('How many elements would you like in the x-direction?\n');
ny = input('How many elements would you like in the y-direction?\n');


check = 1;
while check
    element_type = input('Would you like a truss(1) or CST(0) element?\n');
    if element_type
        
        check = 0;
    elseif element_type == 0
        
        check = 0;
    else
        fprintf('Please enter a valid choice (1 or 0)\n');
    end
    
end
