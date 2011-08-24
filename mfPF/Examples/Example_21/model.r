#########################################################

bolus_1cpt_Vk<-function()
{
	form1<-paste("dose/V*(exp(-k*t))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

immed_lin_null<-function(fconc)
{
  form1<-paste("Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

	formA <- bolus_1cpt_Vk()[[1]]
	formB <- immed_lin_null(formA)[[1]]

form<-c(formA, formB)
