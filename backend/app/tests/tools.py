import logging
import json
from django.urls import reverse
import time
from django.test import TestCase

def logprint(obj, level='error', prettify=True):
    if prettify:
        if isinstance(obj, dict):
            obj = json.dumps(obj, indent=2)
    level = {
        'error': logging.ERROR
    }[level]
    logging.log(level, str(obj))

class AppTestCase(TestCase):
    endpoint = 'endpoint'
    defaults = None

    def run_job(self, job_timeout=20, **params):
        # params['synchronous_job'] = True
        defaults = {} if self.defaults is None else self.defaults
        params_with_defaults = {k: v for (k, v) in defaults.items()}
        params_with_defaults.update(params)
        response = self.client.post(reverse(self.endpoint),
                                    params_with_defaults)
        job_id = response.json()['job_id']
        t0 = time.time()
        while True:
            time.sleep(0.2)
            elapsed = time.time() - t0
            if elapsed > job_timeout:
                raise IOError("Timeout of job at %s" % self.endpoint)
            response = self.client.post(reverse('poll'),
                                        {'job_id': job_id}).json()
            logprint(response)
            if response.get('result', False):
                return response['result']
