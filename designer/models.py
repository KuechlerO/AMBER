from django.db import models



class Alpha_missense(models.Model):
    id = models.IntegerField(primary_key=True)
    protein_variant = models.CharField(max_length=128,blank=True, null=True)
    uniprot_id = models.CharField(max_length=128,blank=True, null=True)
    am_pathogenicity = models.IntegerField()
    am_class = models.CharField(max_length=128,blank=True, null=True)
    class Meta:
        managed = False
        db_table = '"alpha_missense_db"."alpha_missense"'

    # Example
    # uniprot_id      protein_variant am_pathogenicity        am_class
    # A0A024R1R8      M1A             0.4673                  ambiguous
