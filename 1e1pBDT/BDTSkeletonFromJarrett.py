from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from numpy import asarray

################################################################################################################################
################################################################################################################################
## First bit is just looping over the FinalVertexVariable TTrees to fill the training data in
################################################################################################################################

# Signal variable list & label
X0 = []
Y0 = []

# Background variable list & label
X1 = []
Y1 = []

for x in FVV_FOR_SIGNAL_SAMPLE:

	CUTS NEEDED TO DEFINE SIGNAL, e.g. x.is1l1p0pi = True etc

	eventVariables  = [x.Enu_1e1p,x.Electron_Edep,x.Eta,x.PT_1e1p,x.AlphaT_1e1p,x.SphB_1e1p,x.PzEnu_1e1p,x.ChargeNearTrunk*CorrectionFactorPoint(x.Xreco,x.Yreco,x.Zreco),x.Q0_1e1p,x.Q3_1e1p,x.Thetas,x.Phis,x.PTRat_1e1p,x.Proton_TrackLength,x.Lepton_TrackLength,x.Proton_ThetaReco,x.Proton_PhiReco,x.Lepton_ThetaReco,x.Lepton_PhiReco,max(x.MinShrFrac,-1),max(x.MaxShrFrac,-1),x.shower1_smallQ_Y/x.shower1_sumQ_Y ,x.OpenAng]

	X0.append(eventVariables)
	Y0.append(0) #arbitrarily choose 0 as the 'signal' label. Totally unimportant what you pick




for X in FVV_FOR_BACKGROUND:

	eventVariables  = [x.Enu_1e1p,x.Electron_Edep,x.Eta,x.PT_1e1p,x.AlphaT_1e1p,x.SphB_1e1p,x.PzEnu_1e1p,x.ChargeNearTrunk*CorrectionFactorPoint(x.Xreco,x.Yreco,x.Zreco),x.Q0_1e1p,x.Q3_1e1p,x.Thetas,x.Phis,x.PTRat_1e1p,x.Proton_TrackLength,x.Lepton_TrackLength,x.Proton_ThetaReco,x.Proton_PhiReco,x.Lepton_ThetaReco,x.Lepton_PhiReco,max(x.MinShrFrac,-1),max(x.MaxShrFrac,-1),x.shower1_smallQ_Y/x.shower1_sumQ_Y ,x.OpenAng]

	X1.append(eventVariables)



################################################################################################################################
## Loop over as many trees as you want to feed it examples of background
################################################################################################################################
################################################################################################################################


####################################
## Now the actual BDT stuff
####################################


#Class 1 Variables will now have format...
X0  = [ [ event 1 var 1, event 1 var 2, ... event 1 var N]
	[ event 2 var 1, event 2 var 2, ... event 2 var N]
	...
	[ event N var 1, event N var 2, ... event N var N]
      ]

#Class 2 Variables will now have format
X1  = [ [ event 1 var 1, event 1 var 2, ... event 1 var N]
	[ event 2 var 1, event 2 var 2, ... event 2 var N]
	...
	[ event N var 1, event N var 2, ... event N var N]
      ]


# Can generalize this to more classes and do a multi class BDT, or whatever you want.

# Initialize the training & test sample. train_test_split will eat your lists and labels above and split 
# them up for you randomly

seed = 13
test_size = 0.5
x_train, x_test, y_train, y_test = train_test_split(asarray(X0+X1), asarray(Y0+Y1), test_size=test_size, random_state=seed)

# XGBClassifier has a huge number of parameters and there is not honestly a clear cut 'best' set for any 
# particular situation. A large max_depth and large n_estimators will make your training deep and relatively complicated
# Trade off between efficacy and potential for over fitting to simulation. Regularize with gamma. I use deep, large 
# networks for the 1e1p with a fairly large gamma value. You can play with these below. 
# I DON'T KNOW HOW MANY CORES YOU HAVE, DON'T LEAVE IT AT 12 IF YOU DON'T HAVE 12. IT WILL CRASH YOU

MyModel = XGBClassifier(silent=True, 
                       learning_rate=0.1,  
                       objective="binary:logistic",
                       n_estimators=5000,
                       max_depth=18,
                       gamma=1,
                       nthread = 12)

MyModel.fit(x_train, y_train , sample_weight = w_train)


# Depending on how deep and complicated your training is and how many cores you decided to call this might take between 30s and 20 mins.
#Once you havec a trained MyModel, you can get the probability array for an arbitrary MyUnknownEvent = [var1, var2, var3 ...] by calling

MyProbabilities = MyModel.predict_proba(MyUnknownEvent)

#This will return an array of 2 probabilities since you trained on 2 classes. 



# There are some pretty simple adjustments or augmentations I've omitted for simplicity, but if you aren't getting reasonable results we can try those




