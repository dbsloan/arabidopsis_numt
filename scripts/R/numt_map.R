library(diagram)

numt_max = 640560
mito_max = 366924
mito_offset = (1 - mito_max/numt_max)/2

blast_coords = read.csv("numt_map/blast_coords.csv")
bac_coords = read.csv("numt_map/bac_coords.csv")
var_coords = read.csv("numt_map/var_coords.csv")
snp_links = read.csv("numt_map/SNP_links.csv")
feature_coords = read.csv("numt_map/feature_coords.csv")
single_copy = read.csv("numt_map/single_copy.csv")


height1 = 0.55
height2 = 0.85
chrom_nudge = 0.005
chrom_offset1 = -0.09
chrom_offset2 = -0.18

single_copy_offset = 0.015

var_height = 0.025

snp_curve_height = 0.1

repeat_offset = 0.1
repeat_height = 0.03
left1a = 8207/numt_max
right1a = 83859/numt_max
left1b = 183897/numt_max
right1b = 235429/numt_max
left2a = 242658/numt_max
right2a = 318224/numt_max
left2b = 319222/numt_max
right2b = 370754/numt_max
left3a = 377981/numt_max
right3a = 453599/numt_max
left3b = 454597/numt_max
right3b = 505929/numt_max

bac_spacer1 = 0.15
bac_spacer2 = 0.01
bac_height = 0.03
bac_offset = -0.06
bac_title_offset = -0.15
bac_count = 4


legend_offset = 0.03
legend_offset_top = 0.14

tick_length = 0.02
text_offset = 0.02
legend_text_offset = 0.06

feature_height = 0.05
extra_feature_height = 0.035

tandem_center = 335
tandem_length = 135
tandem_offset = 0.09
tandem_height = 0.03
tandem_text_offset = 0.02

missing1_left = 69277
missing1_right = 251962
missing2_left = 327775
missing2_right = 459381
missing_offset = 0.04


plot.new()

for (i in 1:dim(blast_coords)[1]){
  left1 = blast_coords[i, "StartA"]/numt_max
  right1 = blast_coords[i, "EndA"]/numt_max
  left2 = blast_coords[i, "StartB"]/numt_max + mito_offset
  right2 = blast_coords[i,"EndB"]/numt_max + mito_offset
  polygon(c(left1, left2, right2, right1), c(height1, height2, height2, height1), col = adjustcolor("darkolivegreen3", alpha.f=.35), lwd = 0.25)
} 

for (i in 1:dim(single_copy)[1]){
  left = mito_offset + single_copy[i, "Start"]/numt_max
  right = mito_offset + single_copy[i, "End"]/numt_max
  segments(left, height2 + single_copy_offset, right, height2 + single_copy_offset, lwd=3, col=single_copy[i, "Color"])
} 

segments (left1a, height1 - repeat_offset - repeat_height/2, right1b, height1 - repeat_offset - repeat_height/2, col="lightgray")
polygon (c(left1a, left1a, right1a, right1a), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="lightgray", border=NA)
polygon (c(left1b, left1b, right1b, right1b), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="lightgray", border=NA)
#segments (left2a, height1 - repeat_offset - repeat_height/2, right2b, height1 - repeat_offset - repeat_height/2, col="gray")
polygon (c(left2a, left2a, right2b, right2b), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="darkgray", border=NA)
#polygon (c(left2b, left2b, right2b, right2b), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="darkgray", border=NA)
#segments (left3a, height1 - repeat_offset - repeat_height/2, right3b, height1 - repeat_offset - repeat_height/2, col="gray")
polygon (c(left3a, left3a, right3b, right3b), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="black", border=NA)
#polygon (c(left3b, left3b, right3b, right3b), c(height1 - repeat_offset, height1 - repeat_offset - repeat_height, height1 - repeat_offset - repeat_height, height1 - repeat_offset), col="black", border=NA)

text (bac_title_offset, height1 - repeat_offset - repeat_height/2, "3-Copy Reps", adj = 0, cex = 0.75, font=2)


for (i in 1:dim(bac_coords)[1]){
  left = bac_coords[i, "Start"]/numt_max
  right = bac_coords[i, "End"]/numt_max
  bac_height1 = height1 - bac_spacer1 - (bac_coords[i, "Level"] - 1)*(bac_spacer2+bac_height)
  bac_height2 = height1 - bac_spacer1 - (bac_coords[i, "Level"] - 1)*(bac_spacer2+bac_height) - bac_height
  polygon(c(left, left, right, right), c(bac_height1, bac_height2, bac_height2, bac_height1), col = adjustcolor(bac_coords[i, "Color"], alpha.f=bac_coords[i, "Alpha"]), border=NA)
  text (bac_offset, (bac_height1+bac_height2)/2, bac_coords[i, "BAC"], adj = 0, cex = 0.5)
} 
text (bac_title_offset, height1 - bac_spacer1 - ((bac_count - 1)/2)*(bac_spacer2+bac_height) , "BACs", adj = 0, cex = 0.75, font=2)

for (i in 1:dim(snp_links)[1]){
  if (snp_links[i, "Exclude"] == 0){
    SNP1 = snp_links[i, "SNP1"]/numt_max;
    SNP2 = snp_links[i, "SNP2"]/numt_max;
    curvedarrow(c(SNP1, height1 - var_height), c(SNP2, height1 - var_height), arr.type="none", curve = 0.05/(SNP2-SNP1), lwd=0.5, lcol=adjustcolor("gray", alpha.f=0.6))
  }
} 

for (i in 1:dim(var_coords)[1]){
  pos = var_coords[i, "Position"]/numt_max
  segments(pos, height1, pos, height1 - var_height, lwd=0.5, col=adjustcolor(var_coords[i, "Color"], alpha.f = 0.6))
} 


for (i in 1:dim(feature_coords)[1]){
  left = mito_offset + feature_coords[i, "Start"]/numt_max
  right = mito_offset + feature_coords[i, "End"]/numt_max
  polygon(c(left, left, right, right), c(height2 + chrom_nudge, height2 + chrom_nudge + feature_height + extra_feature_height*feature_coords[i, "Extra"], height2 + chrom_nudge + feature_height + extra_feature_height*feature_coords[i, "Extra"], height2 + chrom_nudge), col = "gray85", lwd = 0.5)
  text ((left+right)/2, height2 + chrom_nudge + feature_height + extra_feature_height*feature_coords[i, "Extra"] + text_offset, feature_coords[i, "Feature"], cex = 0.5)
} 

segments(0, height1 - chrom_nudge, 1, height1 - chrom_nudge, lwd=3)
segments(mito_offset, height2 + chrom_nudge, mito_offset + mito_max/numt_max, height2 + chrom_nudge, lwd=3)

text (chrom_offset1, height1 - chrom_nudge, "numt", adj=0, cex=0.75, font=2)
text (mito_offset + chrom_offset2, height2 + chrom_nudge, "Mitogenome", adj=0, cex=0.75, font=2)


tandem_vertical = height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset
tandem_left_pos = 1000*(tandem_center-tandem_length)/numt_max
tandem_center_pos = 1000*tandem_center/numt_max
tandem_right_pos = 1000*(tandem_center+tandem_length)/numt_max

polygon(c(tandem_left_pos, tandem_left_pos, tandem_center_pos, tandem_center_pos), c(tandem_vertical, tandem_vertical + tandem_height, tandem_vertical + tandem_height, tandem_vertical), col = "white", lwd = 1)
polygon(c(tandem_center_pos, tandem_center_pos, tandem_right_pos, tandem_right_pos), c(tandem_vertical, tandem_vertical + tandem_height, tandem_vertical + tandem_height, tandem_vertical), col = "white", lwd = 1)
text (bac_title_offset, tandem_vertical+tandem_text_offset, "Tandem Dup", adj = 0, cex = 0.75, font=2)

missing_vertical = tandem_vertical - missing_offset
missing1_left_pos = missing1_left/numt_max
missing1_right_pos = missing1_right/numt_max
missing2_left_pos = missing2_left/numt_max
missing2_right_pos = missing2_right/numt_max

segments (missing1_left_pos, missing_vertical, missing1_right_pos, missing_vertical, lwd=2)
segments (missing2_left_pos, missing_vertical, missing2_right_pos, missing_vertical, lwd=2)
text (bac_title_offset, missing_vertical, "Omitted Seq", adj = 0, cex = 0.75, font=2)

segments(0, height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset, 1, height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset, lwd=1)
segments(mito_offset, height2 + legend_offset_top, mito_offset + mito_max/numt_max, height2 + legend_offset_top, lwd=1)

for (i in seq(0, numt_max/1000, by = 50)){
  segments (i*1000/numt_max, height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset, i*1000/numt_max, height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset - tick_length, lwd = 1)
  text (i*1000/numt_max, height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset - tick_length - text_offset, cex=0.5, i)    
}

text (0.5,  height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset - tick_length - text_offset - legend_text_offset, cex=0.75, font=2, "Position (kb)")
text (0,  height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset - tick_length - text_offset - legend_text_offset, cex=0.75, font=1, adj=0, "Telomere End")
text (1,  height1 - bac_spacer1 - (bac_count - 1)*(bac_spacer2+bac_height) - tandem_offset - missing_offset - legend_offset - tick_length - text_offset - legend_text_offset, cex=0.75, font=1, adj=1, "Centromere End")



for (i in seq(0, mito_max/1000, by = 50)){
  segments (mito_offset + i*1000/numt_max, height2 + legend_offset_top, mito_offset + i*1000/numt_max, height2 + legend_offset_top + tick_length, lwd = 1)
  text (mito_offset + i*1000/numt_max, height2 + legend_offset_top + tick_length + text_offset, cex=0.5, i)    
}  
