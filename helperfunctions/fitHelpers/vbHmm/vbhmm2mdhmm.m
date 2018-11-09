%VBHMM2MDHMM Add minimum duration constraint to standard HMM
%
%        MODEL = VBHMM2MDHMM(MODEL,MD)
%
% Convert a VBHMM to a Minimum Duration VBHMM, with a minimum
% duration of MD. 


function model = hmm2mdhmm(model,md)
if nargin<2
	md = [];
end

% Check:
if isempty(md) || (md<=1)
	warning('moghmm:minDurationOfOne',...
        'Minimum duration is 1, now using the standard HMM.');
	model = hmm;
	model.md = 1;
	return;
end

% initialize:
Q = length(hmm.prior);

% prior is easy:
model.prior = zeros(Q*md,1);
model.prior(1:md:Q*md,1) = hmm.prior;

% transition is harder:
% first create a 'diagonal' matrix where the diagonal is one element
% off:
z1 = zeros(Q*md,1);
z2 = z1; z2(2)=1;
trans = toeplitz(z1,z2);
% then copy the elements on the right spot:
for i=1:Q
	trans(i*md,i*md) = hmm.trans(i,i);
	for j=i+1:Q
		trans(i*md,(j-1)*md+1) = hmm.trans(i,j);
		trans(j*md,(i-1)*md+1) = hmm.trans(j,i);
	end
end
model.trans = sparse(trans);

% % all the pdf's do not change;
% model.pdf = hmm.pdf;

% finally add the fact that we are working with minimum duration:
model.md = md;


