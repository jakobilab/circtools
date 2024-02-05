circrna_length = int(entry[3]) - int(entry[2])

exon1_length = len(exon_cache[circular_rna_id][1])
exon2_length = len(exon_cache[circular_rna_id][2])

exon2_colour = "#ffac68"

if exon2_length == 0:
    exon1_length = int(len(exon_cache[circular_rna_id][1])/2)+1
    exon2_length = int(len(exon_cache[circular_rna_id][1])/2)
    exon2_colour = "#ff6877"

index = int(entry[5])

gdd = GenomeDiagram.Diagram('circRNA probe diagram')
gdt_features = gdd.new_track(1, greytrack=True, name="", )
gds_features = gdt_features.new_set()
 
# adding first exon
feature = SeqFeature(FeatureLocation(0, exon1_length))
feature.location.strand = +1
gds_features.add_feature(feature, name="Exon 1", label=False, color="#ff6877", label_size=22)
# adding second exon
feature = SeqFeature(FeatureLocation(circrna_length - exon2_length, circrna_length))
feature.location.strand = +1
gds_features.add_feature(feature, name="Exon 2", label=False, color="#ffac68", label_size=22)

# # adding the product of 50bp
# feature = SeqFeature(FeatureLocation(0, 25))
# feature.location.strand = -1
# gds_features.add_feature(feature, name="Product", label=False, color="blue")
# feature = SeqFeature(FeatureLocation(circrna_length-25, circrna_length))
# feature.location.strand = -1
# gds_features.add_feature(feature, name="Product: " + str(product_size) + "bp", label=False, color="blue", label_size=22, label_position="middle")

# adding the individual arm
flag = 0            # flag to keep track if start of RBD3 is before BSJ
rbd5_start = circrna_length-25+index 
rbd5_end = circrna_length-25+index+20
if (rbd5_end > circrna_length):
feature = SeqFeature(FeatureLocation(rbd5_start, circrna_length)) 
gds_features.add_feature(feature, name="RBD5", label=False, color="red", label_size=22)
feature = SeqFeature(FeatureLocation(0, rbd5_end - circrna_length))
gds_features.add_feature(feature, name="RBD5", label=False, color="red", label_size=22)
rbd5_end = rbd5_end - circrna_length

else:
feature = SeqFeature(FeatureLocation(rbd5_start, rbd5_end))
gds_features.add_feature(feature, name="RBD5", label=False, color="red", label_size=22)
flag = 1
rbd3_start = rbd5_end + 1

if (flag == 1):
feature = SeqFeature(FeatureLocation(rbd3_start, circrna_length)) 
gds_features.add_feature(feature, name="RBD3", label=False, color="green", label_size=22)
feature = SeqFeature(FeatureLocation(0, 20 - (circrna_length - rbd3_start)))
gds_features.add_feature(feature, name="RBD3", label=False, color="green", label_size=22)
else:
rbd3_end = rbd3_start + 20
feature = SeqFeature(FeatureLocation(rbd3_start, rbd3_end))
gds_features.add_feature(feature, name="RBD3", label=False, color="green", label_size=22)

feature = SeqFeature(FeatureLocation(0, 1))
gds_features.add_feature(feature, name="BSJ", label=True, color="white", label_size=22)
gdd.draw(format='circular', pagesize=(600, 600), circle_core=0.25, track_size=0.2, tracklines=0, x=0.00, y=0.00, start=0, end=circrna_length-1)
gdd.write("test.svg", "SVG")