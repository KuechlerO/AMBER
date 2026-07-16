import os
import sys

from django.apps import AppConfig


class DesignerConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'designer'

    def ready(self):
        # runserver autoreloader starts apps.ready() in parent + child; only the child serves HTTP.
        if 'runserver' in sys.argv and os.environ.get('RUN_MAIN') != 'true':
            return
        try:
            from designer.screen_plot.gene_lookup import is_screen_library_data_readable
            if is_screen_library_data_readable():
                from designer.screen_plot.data_loader import get_screen_data_store
                get_screen_data_store().start_background_load()
        except Exception:
            pass
