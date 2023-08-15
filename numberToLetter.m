function letter = numberToLetter(num)
% Converts a positive integer to a letter sequence using the English alphabet.
% 
% Inputs:
%   num - Positive integer to be converted.
%
% Outputs:
%   letter - Letter sequence representing the input number.
%
% Description:
% The numberToLetter function converts a given positive integer into a letter
% sequence using the English alphabet. It represents each digit in base 26
% (corresponding to the 26 letters in the alphabet) and recursively breaks
% down the number to generate the letter sequence.
%
% Examples:
%   1. numberToLetter(1) returns 'a'
%   2. numberToLetter(26) returns 'z'
%   3. numberToLetter(27) returns 'aa'
%   4. numberToLetter(52) returns 'az'
%   5. numberToLetter(53) returns 'ba'
%   6. numberToLetter(701) returns 'zy'
%   7. numberToLetter(702) returns 'zz'
%   8. numberToLetter(703) returns 'aaa'
    
    base = 26;  % Number of letters in the alphabet
    letters = 'abcdefghijklmnopqrstuvwxyz';
    
    if num <= base
        letter = letters(num);
    else
        quotient = floor((num - 1) / base);
        remainder = mod(num - 1, base) + 1;
        letter = [numberToLetter(quotient), letters(remainder)];
    end
end
