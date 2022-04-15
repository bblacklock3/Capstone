function out = isKicking(hipAng, hipVel)
rom = max(hipAng) - min(hipAng);
normRom = rom / 136; % 136 degrees is average range of motion for hip flexion in infants under 2: https://musculoskeletalkey.com/pediatric-range-of-motion/
normVel = 136/.41; % because an average kick flexion phase is .41 seconds for a 3 month old: https://journals.sagepub.com/doi/full/10.1177/2055668317717461
out = normRom + normVel > 1;
end