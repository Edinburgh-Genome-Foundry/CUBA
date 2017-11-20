<template lang='pug'>
div
  input(type='text', size='25', value='', v-model='honeypot', :style="{display: 'none'}")
  el-button.center(v-if='submitButton && (!polling.inProgress)' @click='submit()') {{submitButtonText }}
  .polling(v-if='polling.inProgress')
    spinner(color="#6da5ff" size="12px")
    .polling-message {{polling.message}}
    el-steps(:space='100', :active='progressStage')
      el-step(icon='upload')
      el-step(icon='setting')
      el-step.last-step(icon='message')

  //- el-alert(v-if='result.error', :title="result.error", type="warning")
  el-alert(v-if='requestError', :title="requestError", type="error", :closable="false")
  .results(v-if='!polling.inProgress')
    downloadbutton(:filedata='result.file' v-if='result.file')
    .results-summary(v-if='result.preview', v-html="result.preview.html")
</template>

<script>
import downloadbutton from './DownloadButton'
import spinner from 'vue-spinner/src/PulseLoader'

export default {
  props: {
    submitButton: {default: true},
    submitButtonText: {default: 'Submit'},
    value: {default: () => ({})},
    backendUrl: {default: ''},
    backendIP: {default: 'auto'},
    form: {default: () => ({})},
    validate_form: {default: () => () => ([])}
  },
  data () {
    console.log(this.backendIP === 'auto' ? this.computeBackendIP() : this.backendIP)
    return {
      honeypot: '',
      backendRoot: this.backendIP === 'auto' ? this.computeBackendIP() : this.backendIP,
      polling: {
        inProgress: false,
        status: 'queued',
        message: '',
        data: ''
      },
      result: {},
      requestError: ''
    }
  },
  components: {
    downloadbutton,
    spinner
  },
  watch: {
    form: {
      handler (val) {
        console.log(val)
        this.$emit('input', val)
      },
      deep: true
    },
    value: {
      handler (val) {
        this.form = val
      },
      deep: true
    }
  },
  computed: {
    fulldata () {
      return {
        data: this.form,
        honeypot: this.honeypot
      }
    },
    progressStage () {
      return ['queued', 'started', 'finished'].indexOf(this.polling.status) + 1
    }
  },
  methods: {
    submit () {
      var errors = this.validate_form()
      if (errors.length) {
        this.requestError = 'Invalid form: ' + errors.join('   ')
        return false
      }
      this.$http.post(
        this.backendRoot + this.backendUrl,
        this.form
      ).then(function (response) {
        // SUCCESS
        console.log('job-start response', response)
        this.startPolling(response.body.job_id)
      }, function (response) {
        // FAILURE OF THE JOB STARTING
        console.log('job-start error response', response)
        this.requestError = 'Failed with status ' + response.status
        if (response.status === 400) {
          this.requestError = 'Bad Request: ' + window.JSON.stringify(response.body)
        }
      })
    },
    startPolling (jobId) {
      this.polling.inProgress = true
      var jobPoller = setInterval(function () {
        this.$http.post(
          this.backendRoot + 'poll',
          {job_id: jobId}
        ).then(function (pollResponse) {
          // console.log('poll response', pollResponse)
          this.requestError = pollResponse.error
          var data = pollResponse.body
          this.polling.status = data.status
          if (data.success === false) {
            this.requestError = data.error
            this.polling.inProgress = false
            clearInterval(jobPoller)
          } else if (data.status === 'queued') {
            this.polling.message = 'Job pending...'
          } else if (data.status === 'started') {
            this.polling.message = data.progress.message
            this.polling.data = data.progress.data
          } else if (data.status === 'finished') {
            this.polling.message = 'Finished ! sending now'
            this.result = data.result
            this.polling.inProgress = false
            clearInterval(jobPoller)
          }
          this.polling.message = data.progress.message
        }, function (pollResponse) {
          // FAILURE OF THE POLLING
          this.polling.inProgress = false
          this.requestError = 'Failed with status ' + pollResponse.status
          clearInterval(jobPoller)
        })
      }.bind(this), 100)
    },
    computeBackendIP () {
      var location = window.location.origin
      if (location[location.length - 5] === ':') {
        location = location.slice(0, location.length - 5)
      }
      return location + ':8081/'
    }
  }
}
</script>

<style scoped>

.submit-button {
  margin-top: 20px;
  margin-bottom: 60px;
  width: 150px;
}

.el-alert {
  margin: 15px auto 15px;
  width: auto;
  font-size: 14px;
  max-width: 700px;
  padding-right: 40px
}

.polling {
  text-align:center;
  margin-top: 50px
}

.last-step {
  width: 40px !important;
}

</style>
