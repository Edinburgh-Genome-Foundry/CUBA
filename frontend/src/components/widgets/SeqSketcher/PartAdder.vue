<template lang='pug'>
.sketcher-part.adder
  p.plus.construct-hover-only +
  .symbol.none(@click="$emit('click')")
</template>

<script>
import part from './Part'

export default {
  name: 'part-adder',
  mixins: [part]
}
</script>
<style lang='scss'>
.sketcher-part.adder {
  min-width: 1.5em;
  text-align: center;
  p.plus {
    min-width: 1em;
    font-size: 2.5em;
    top: 1.3em;
    color: #ccc;
    cursor: pointer;
    padding: 0;
    margin-left: 0;
    margin-right: 0;
    text-align: center;

  }
}
</style>
