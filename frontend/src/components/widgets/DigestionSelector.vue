<template lang="pug">
el-select(v-model='selectedEnzymes', placeholder='Enzymes in this digestion', filterable, :multiple='true')
  el-option(v-for='enzyme in enzymes.list', :label='enzyme', :value='enzyme', :key='enzyme')
</template>

<script>
import enzymes from '../../assets/enzymes'

export default {
  props: {
    value: {
      default: () => ([])
    }
  },
  data () {
    return {
      enzymes: enzymes,
      selectedEnzymes: this.value
    }
  },
  watch: {
    value: {
      deep: true,
      handler (val) {
        this.selectedEnzymes = val
      }
    },
    selectedEnzymes: {
      deep: true,
      handler (val) {
        this.$emit('input', val)
      }
    }
  }
}
</script>

<style scoped>
.el-tag {
  height: 35px !important;
}

.el-select {
  width: 100%;
}
</style>
