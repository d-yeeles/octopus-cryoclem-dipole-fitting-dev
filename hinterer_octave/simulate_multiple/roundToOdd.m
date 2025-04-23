function rounded = roundToOdd(f)
    rounded = round(f); % Round to the nearest integer
    if mod(rounded, 2) == 0 % Check if it's even
        rounded = rounded - 1 * sign(f); % Adjust to the nearest odd integer
    end
end