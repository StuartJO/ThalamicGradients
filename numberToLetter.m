function letter = numberToLetter(num)
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