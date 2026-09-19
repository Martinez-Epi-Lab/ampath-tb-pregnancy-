# Manuscript figures from aggregate analysis results.
figure_palette <- function() c(reference='#E5E9EC',pregnancy='#28617A',postpartum='#80B8C5',excluded='#A7ADB3',cohort='#235B79',sccs='#9C514A',ink='#243440')

save_manuscript_figure <- function(directory,name,width,height,draw) {
  paths<-setNames(file.path(directory,paste0(name,c('.pdf','.svg','.png'))),c('pdf','svg','png'))
  for(type in names(paths)){
    if(type=='pdf')grDevices::cairo_pdf(paths[[type]],width=width,height=height,family='Helvetica',onefile=FALSE)
    else if(type=='svg')grDevices::svg(paths[[type]],width=width,height=height,family='sans',onefile=FALSE)
    else grDevices::png(paths[[type]],width=width,height=height,units='in',res=300,type='cairo',family='sans')
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  paths
}

figure_study_design <- function() {
  pal<-figure_palette();par(mar=c(3.2,3,2.2,1),mgp=c(2,.65,0),family='sans',fg=pal['ink'],col.axis=pal['ink'])
  plot.new();plot.window(xlim=c(-130,1140),ylim=c(0,10),xaxs='i',yaxs='i')
  title('Observation and exposure periods',adj=0,cex.main=1.25)
  bars<-function(y,stop,full){
    rect(0,y-.37,stop,y+.37,col=pal['reference'],border=NA)
    rect(0,y-.37,80,y+.37,col=pal['pregnancy'],border=NA)
    rect(81,y-.37,min(150,stop),y+.37,col=pal['postpartum'],border=NA)
    rect(151,y-.37,min(360,stop),y+.37,col=pal['pregnancy'],border=NA)
    if(stop>=361)rect(361,y-.37,min(540,stop),y+.37,col=pal['postpartum'],border=NA)
    if(full){rect(650,y-.37,970,y+.37,col='#F6F7F8',border=pal['excluded'],density=13,angle=45);segments(790,y-.37,790,y+.37,col=pal['excluded'])}
    rect(-100,y-.37,0,y+.37,col='#FFFFFF',border=pal['pregnancy'],lty=2)
    segments(0,y-.62,0,y+.62,col=pal['ink'],lty=3)
    points(450,y,pch=18,col='#A13D38',cex=1.5)
  }
  text(-125,8.7,'A  Cohort',adj=0,font=2,cex=1.12);bars(7.6,450,FALSE)
  segments(450,7.6,1100,7.6,lty=3,col='#939CA2')
  text(770,8.05,'Follow-up stops at first TB event',cex=.92)
  text(40,7.6,'1',col='white',font=2);text(255,7.6,'2',col='white',font=2)
  text(450,6.75,'TB event',col='#A13D38',cex=.9)
  text(-125,5.6,'B  Self-controlled case series',adj=0,font=2,cex=1.12);bars(4.5,1100,TRUE)
  text(40,4.5,'1',col='white',font=2);text(255,4.5,'2',col='white',font=2)
  text(810,5.35,'Third pregnancy and postpartum excluded',cex=.86)
  text(450,3.68,'TB event',col='#A13D38',cex=.9)
  text(1028,3.68,'Later reference\ntime retained',cex=.83)
  arrows(115,6.75,115,7.12,length=.07,col=pal['ink']);text(115,6.4,'Postpartum ends at next pregnancy',cex=.81)
  text(-110,2.8,'Pregnancy 1 began\nbefore entry',adj=0,cex=.82)
  segments(0,2.1,1100,2.1,col=pal['ink']);for(t in c(0,180,360,540,720,900,1080)){segments(t,2.1,t,1.98);text(t,1.65,t,cex=.82)}
  text(550,1.13,'Illustrative days since study entry',cex=.94)
  legend('bottom',inset=c(0,-.005),horiz=TRUE,bty='n',fill=pal[c('reference','pregnancy','postpartum','excluded')],border=NA,
    legend=c('Reference time','Pregnancy','Postpartum','Excluded time'),cex=.84,x.intersp=.65)
}

figure_primary_forest <- function(coefficients) {
  pal<-figure_palette();x<-coefficients[model_id %in% c('cohort_primary','sccs_primary') & feature %in% c('pregnancy','postpartum_days_1_180')]
  x[,order_row:=match(paste(model_id,feature),c('cohort_primary pregnancy','cohort_primary postpartum_days_1_180','sccs_primary pregnancy','sccs_primary postpartum_days_1_180'))];data.table::setorder(x,order_row)
  par(mar=c(4,1,2.6,1),family='sans',fg=pal['ink']);plot.new();plot.window(xlim=c(0,1),ylim=c(0,5.8),xaxs='i',yaxs='i')
  title('Primary associations with tuberculosis',adj=0,cex.main=1.2)
  map<-function(v).38+(log(v)-log(.7))/(log(2.5)-log(.7))*.37
  segments(map(1),.6,map(1),4.9,lty=2,col='#9CA5AB')
  text(.01,5.28,'Study design / exposure',adj=0,font=2,cex=.96);text(.99,5.28,'Adjusted IRR (95% CI)',adj=1,font=2,cex=.96)
  yy<-c(4.25,3.45,2,1.2)
  text(.01,4.77,'Cohort',adj=0,font=2,col=pal['cohort']);text(.01,2.52,'SCCS',adj=0,font=2,col=pal['sccs'])
  for(i in seq_len(nrow(x))){col<-if(i<=2)pal['cohort']else pal['sccs'];text(.03,yy[i],if(x$feature[i]=='pregnancy')'Pregnancy'else'Postpartum days 1-180',adj=0,cex=.92)
    segments(map(x$ci_lower[i]),yy[i],map(x$ci_upper[i]),yy[i],lwd=2,col=col);points(map(x$irr[i]),yy[i],pch=15,cex=1.05,col=col)
    text(.99,yy[i],sprintf('%.2f (%.2f-%.2f)',x$irr[i],x$ci_lower[i],x$ci_upper[i]),adj=1,cex=.95)
  }
  ticks<-c(.75,1,1.5,2,2.5);axis(1,at=map(ticks),labels=ticks,pos=.48,tck=-.018,cex.axis=.85)
  text(.565,-.5,'Incidence rate ratio (log scale)',xpd=NA,cex=.93)
  mtext('Age- and calendar-adjusted; first two documented pregnancies.',side=1,line=2.8,adj=0,cex=.8)
}

figure_participant_flow <- function(flow) {
  pal<-figure_palette();n<-setNames(flow$count,flow$key);f<-function(key)format(n[[key]],big.mark=',',trim=TRUE,scientific=FALSE)
  par(mar=c(.3,.3,1.7,.3),family='sans',fg=pal['ink']);plot.new();plot.window(xlim=c(0,10),ylim=c(0,10),xaxs='i',yaxs='i');title('Participant and event selection',adj=0,cex.main=1.2)
  box<-function(x,y,w,h,label,main=TRUE){rect(x-w/2,y-h/2,x+w/2,y+h/2,col=if(main)'#EDF3F6'else'#F7F7F7',border=if(main)pal['cohort']else'#B1B8BD',lwd=.9);text(x,y,label,cex=if(main).99 else .88)}
  arrow<-function(x1,y1,x2,y2)arrows(x1,y1,x2,y2,length=.08,col='#667780')
  center<-3.4;side<-8.1;ys<-c(9.1,7.7,6.3,4.9,3.5,2.1)
  labels<-c(paste0('Eligible observation window\n',f('observation_all'),' people'),paste0('Female-coded participants\n',f('observation_female'),' women'),paste0('At least one eligible pregnancy\n',f('pregnancy_eligible'),' women'),paste0('Incident-TB risk cohort\n',f('final_cohort'),' women'),paste0('First incident TB events\n',f('incident_events')),paste0('Primary analyzed TB events / SCCS cases\n',f('analyzed_events')))
  for(i in seq_along(ys)){box(center,ys[i],5.6,.83,labels[i]);if(i<length(ys))arrow(center,ys[i]-.415,center,ys[i+1]+.415)}
  box(side,7.7,3.25,.86,paste0(f('observation_men'),' male-coded\nparticipants'),FALSE);arrow(center,8.42,side,8.42);arrow(side,8.42,side,8.15)
  without<-n[['observation_female']]-n[['pregnancy_eligible']]
  box(side,6.3,3.25,.86,paste0(format(without,big.mark=',',trim=TRUE),' without\nan eligible pregnancy'),FALSE);arrow(center,7,side,7);arrow(side,7,side,6.75)
  box(side,4.9,3.25,1.05,paste0('Prior TB exclusions\n',f('TB_prior_treatment_excluded'),' recorded prior treatment\n',f('TB_broad_evidence_excluded'),' earlier TB evidence'),FALSE);arrow(center,5.6,side,5.6);arrow(side,5.6,side,5.45)
  box(side,3.5,3.25,.86,paste0(format(n[['final_cohort']]-n[['incident_events']],big.mark=',',trim=TRUE),' without\nan incident TB event'),FALSE);arrow(center,4.2,side,4.2);arrow(side,4.2,side,3.95)
  box(side,2.1,3.25,1.02,paste0(f('excluded_later_pregnancy_events'),' events in third-or-later\npregnancy or postpartum\ntime excluded'),FALSE);arrow(center,2.8,side,2.8);arrow(side,2.8,side,2.65)
  text(5,.9,'Pregnancy eligibility and prior-TB exclusions are applied sequentially.',cex=.83)
  text(5,.52,'Primary analyses exclude third-or-later pregnancy and its postpartum time.',cex=.83)
}

figure_postpartum_spline <- function(curves,tests) {
  pal<-figure_palette();par(mfrow=c(1,2),mar=c(4.6,4.4,3,1),oma=c(1,0,0,0),mgp=c(2.5,.7,0),family='sans',fg=pal['ink'])
  for(design in c('cohort','sccs')){
    selected_design<-design
    z<-curves[design==selected_design];data.table::setorder(z,postpartum_day);color<-pal[design]
    plot(z$postpartum_day,z$irr,type='n',log='y',ylim=c(.65,3.2),xlim=c(1,180),axes=FALSE,xlab='Days after pregnancy end',ylab='Adjusted incidence rate ratio',main=if(design=='cohort')'A  Cohort'else'B  SCCS',cex.main=1.1)
    polygon(c(z$postpartum_day,rev(z$postpartum_day)),c(z$ci_lower,rev(z$ci_upper)),col=grDevices::adjustcolor(color,alpha.f=.17),border=NA)
    abline(h=1,lty=2,col='#87939B');lines(z$postpartum_day,z$irr,col=color,lwd=2)
    axis(1,at=c(1,30,60,90,120,150,180));axis(2,at=c(.75,1,1.5,2,3),labels=c('0.75','1.0','1.5','2.0','3.0'),las=1);box(bty='l')
    tt<-tests[design==selected_design];p<-setNames(tt$p_value,tt$test)
    legend('topright',legend=c(sprintf('Overall curve: P = %.3f',p['overall_postpartum_curve']),sprintf('Time variation: P = %.3f',p['variation_over_postpartum_time']),sprintf('Nonlinearity: P = %.3f',p['nonlinearity'])),bty='n',cex=.81,inset=.01)
    rug(c(10,64,117,171),side=1,col=color,ticksize=.025)
  }
  mtext('Shading: 95% confidence interval. Restricted cubic spline knots: days 10, 64, 117 and 171.',side=1,outer=TRUE,cex=.82,line=0)
}
