% List of datasets available for recon.
% Dataset 1: Ok, but smudged periods
% Dataset 2: Ok
% Dataset 3: Not ok, the cable was taut
% Dataset 4: Ok, but  a bit of intial strangeness
% Dataset 5: Ok
% Dataset 6: Great
% Dataset 7: Ok
% Dataset 8: Ok, a lot of cardiac signal

function [S stemdir] = inputdatasets();

stemdir = '/datadrive/OCM/PET/';

% Patient #1, imaged on 2023-06-26
S(01) = struct('dir', [stemdir 'OCM2023-06-26/'], 'fname', 'OCM2023-06-26-data001', 'tag', 'PAT01');

% Patient #2, imaged on 2023-09-01
S(02) = struct('dir', [stemdir 'OCM2023-09-01-A/'], 'fname', 'OCM2023-09-01-data001', 'tag', 'PAT02');

% Patient #3, imaged on 2023-09-01 as well
S(03) = struct('dir', [stemdir 'OCM2023-09-01-B/'], 'fname', 'OCM2023-09-01-data003', 'tag', 'PAT03');

% Patient #4, imaged on 2024-01-16
S(04) = struct('dir', [stemdir 'OCM2024-01-16-A/'], 'fname', 'OCM2024-01-16-data002', 'tag', 'PAT04');

% Patient #5, imaged on 2024-01-16 as well
S(05) = struct('dir', [stemdir 'OCM2024-01-16-B/'], 'fname', 'OCM2024-01-16-data003', 'tag', 'PAT05');

% Patient #6, imaged on 2024-01-16 as well
S(06) = struct('dir', [stemdir 'OCM2024-01-16-C/'], 'fname', 'OCM2024-01-16-data004', 'tag', 'PAT06');

% Patient #7, imaged on 2024-01-16 as well
S(07) = struct('dir', [stemdir 'OCM2024-01-16-D/'], 'fname', 'OCM2024-01-16-data005', 'tag', 'PAT07');

% Patient #8, imaged on 2024-01-18
S(08) = struct('dir', [stemdir 'OCM2024-01-18/'], 'fname', 'OCM2024-01-18-data001', 'tag', 'PAT08');

